#include <iostream>
#include <numbers>

#include <Injector/Grids/Defines.h>
#include <Injector/Grids/GridsFactory.hpp>
#include <Injector/Properties/LogsFactory.hpp>
#include <Injector/Model/Well/Well.hpp>
#include <Injector/Model/Phases/FluidFactory.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


#include "includes/transfer_to_eigen.hpp"
#include "includes/make_steps.hpp"

using namespace Catch;
using namespace Catch::Matchers;

using namespace std;

using namespace GPN;
using namespace GPN::Grids;
using namespace GPN::Logs;
using namespace GPN::Phases;
using namespace GPN::CoordinateTypes;

std::vector<RealType> grid_stencils{0.0, 1.0, 3.0, 7.0, 10.0};
std::vector<RealType> grid_thickness{make_steps(grid_stencils)};
std::vector<RealType> permeability_stencils(grid_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_permeable_stencils(grid_stencils.size() - 1ull, 1.0);
std::vector<RealType> is_perforated_stencils{is_permeable_stencils};

const RealType rate{1.0};
const RealType pressure{1.0 / (2.0 * std::numbers::pi * permeability_stencils.back() * grid_stencils.back())};

const RealType rMax{std::numbers::e}; // m
/*well*/
const RealType sandface_radius{1.0}; // m
const RealType column_radius{0.5}; // m
const RealType tube_radius{0.1};     // m

const RealType tol{1e-12};

TEST_CASE("Well_Test")
{
    is_perforated_stencils[0ll] = 0.0;

    const auto z_grid{
        Grids::Factory::create_axes<CoordinateTypes::Z>(
            grid_stencils)};

    const auto is_permeable{
        Logs::IsPermeableFactory::create(
            is_permeable_stencils,
            z_grid)};

    const auto is_perforated{
        Logs::IsPerforatedFactory::create(
            is_perforated_stencils,
            is_permeable_stencils,
            z_grid)};

    const auto permeability{
        Logs::PermeabilityFactory::create(
            permeability_stencils,
            is_permeable_stencils,
            z_grid)};

    const auto water{FluidFactory::create_water(
        Viscosity{1.0},
        Density{1.0},
        SpecificHeatCapacity{1.0},
        GPN::HeatConductivity{1.0})};

    WellHoles well_holes{
        TubeInnerRadius{tube_radius}, 
        ColumnOuterRadius{column_radius}, 
        SandfaceRadius{sandface_radius}};

    cout << "thickness profile:\n"
         << transfer_to_eigen(grid_thickness).transpose();

    struct Record
    {
        const RealType rate, pressure;
    } history_record_q{rate, std::numeric_limits<double>::quiet_NaN()},
        history_record_p{std::numeric_limits<double>::quiet_NaN(), pressure};

    const Well_KH well_q{
        water, is_permeable, is_perforated, permeability, well_holes, rMax};

    const auto rfp_q = RFPFactory::create_from_container<Logs::RFP>(well_q.get_RFP(history_record_q), is_permeable);
    const auto wfp_q = WFPFactory::create_from_container<Logs::WFP>(well_q.get_WFP(history_record_q), is_perforated);

    const Well_KH well_p{
        water, is_permeable, is_perforated, permeability, well_holes, rMax};

    const auto rfp_p = RFPFactory::create_from_container<Logs::RFP>(well_p.get_RFP(history_record_p), is_permeable);
    const auto wfp_p = WFPFactory::create_from_container<Logs::WFP>(well_p.get_WFP(history_record_p), is_perforated);

    { // check Well_KH

        { // check well_kh at fixed rate

            cout << "\nwell rate profile:     \n"
                 << wfp_q.log_vals.transpose();
            cout << "\nreservoir rate profile:\n"
                 << rfp_q.log_vals.transpose();

            for (auto i{0ll}; i < rfp_q.size(); ++i)
            {
                CHECK_THAT(rate / (grid_stencils.back() - grid_stencils.front()) *
                               z_grid.dual_steps(i),
                           WithinRel(rfp_q(i), tol));
            }

            CHECK_THAT(rfp_q.log_vals.sum(),
                       WithinRel(rate, tol));
            CHECK_THAT(wfp_q.log_vals.sum(),
                       WithinRel(rate, tol));

            RealType cum_rate{0.0};
            ptrdiff_t i{0ll};
            for (auto l{0ll}; (i < rfp_q.size()) && (l < 2ll); ++i)
            {
                if (is_permeable(i) == 1.0)
                {
                    cum_rate += rfp_q(i);
                    ++l;
                }
            }

            CHECK(wfp_q(i - 1ll) == cum_rate);
        }

        { // check well_kh at fixed pressure
            cout << "\nwell rate profile:     \n"
                 << wfp_p.log_vals.transpose();
            cout << "\nreservoir rate profile:\n"
                 << rfp_p.log_vals.transpose();

            for (auto i{0ll}; i < rfp_p.size(); ++i)
            {
                CHECK_THAT(rate / (grid_stencils.back() - grid_stencils.front()) *
                               z_grid.dual_steps(i),
                           WithinRel(rfp_p(i), tol));
            }

            CHECK(rfp_p.log_vals.sum() == rate);
            CHECK(wfp_p.log_vals.sum() == rate);

            RealType cum_rate{0.0};
            ptrdiff_t i{0ll};
            for (auto l{0ll}; (i < rfp_p.size()) && (l < 2ll); ++i)
            {
                if (is_permeable(i) == 1.0)
                {
                    cum_rate += rfp_p(i);
                    ++l;
                }
            }

            CHECK(wfp_p(i - 1ll) == cum_rate);
        }

        { // compare well_kh at fixed rate vs fixed pressure
            assert(rfp_p.size() == rfp_q.size());
            assert(wfp_p.size() == wfp_q.size());
            assert(rfp_p.size() == wfp_p.size());
            assert(rfp_q.size() == rfp_q.size());
            for (auto i{0ll}; i < rfp_p.size(); ++i)
            {
                CHECK_THAT(rfp_p.log_vals(i), WithinRel(rfp_q.log_vals(i), tol));
                CHECK_THAT(wfp_p.log_vals(i), WithinRel(wfp_q.log_vals(i), tol));
            }
        }
    }

    { // check Well_explicit
        vector<RealType> weights;
        weights.reserve(grid_thickness.size());
        for (auto i{0ull}; i < grid_thickness.size(); ++i)
            weights.push_back(grid_thickness[i] * permeability_stencils[i]);

        const Well_Explicit well_q_exp{
            is_permeable, is_perforated, transfer_to_eigen(weights)};

        const auto rfp_q_exp = RFPFactory::create_from_container<Logs::RFP>(well_q_exp.get_RFP(history_record_q), is_permeable);
        const auto wfp_q_exp = WFPFactory::create_from_container<Logs::WFP>(well_q_exp.get_WFP(history_record_q), is_perforated);

        assert(rfp_p.size() == rfp_q_exp.size());
        assert(wfp_p.size() == wfp_q_exp.size());
        assert(rfp_p.size() == wfp_p.size());
        assert(rfp_q_exp.size() == rfp_q_exp.size());
        for (auto i{0ll}; i < rfp_p.size(); ++i)
        {
            CHECK_THAT(rfp_p.log_vals(i), WithinRel(rfp_q_exp.log_vals(i), tol));
            CHECK_THAT(wfp_p.log_vals(i), WithinRel(wfp_q_exp.log_vals(i), tol));
        }
    }
}
