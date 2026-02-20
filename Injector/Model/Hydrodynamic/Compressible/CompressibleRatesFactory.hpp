#pragma once

#include <cassert>
#include <vector>

#include <Injector/Grids/Defines.h>

#include <Injector/History/InjectorRegimes.hpp>

#include <Injector/Properties/FlowField.hpp>
#include <Injector/Properties/FaceProperties.hpp>
#include <Injector/Properties/JT_FieldFactory.hpp>

#include <Injector/Model/Hydrodynamic/SomeRatesFactory.hpp>

namespace GPN
{
    namespace FaceProperties
    {
        template <
            typename Grid2D_t,
            typename Well_t,
            typename History_t,
            typename Fluid_t,
            typename Hydrodynamics_t>
        struct CompressibleRatesFactory : public SomeRatesFactory<Grid2D_t, Well_t, History_t, Fluid_t, Hydrodynamics_t>
        {
            using Base = SomeRatesFactory<Grid2D_t, Well_t, History_t, Fluid_t, Hydrodynamics_t>;
            using Base::first_size;
            using Base::get_history_record;
            using Base::get_pressure_field;
            using Base::history;
            using Base::second_size;
            using Base::solution;

            CompressibleRatesFactory(
                cptr<Hydrodynamics_t> pressure_field,
                const cptr<Grid2D_t> grid2D_rocks,
                const ptr<Well_t> well,
                const ptr<History_t> history,
                const Fluid_t &fluid)
                : Base{pressure_field, grid2D_rocks, well, history, fluid},
                  //  face_mobility{pressure_field->face_mobility},
                  mobility_factor{
                      set_mobility_factor(
                          pressure_field->face_mobility.face_vals_axes2,
                          grid2D_rocks->first_coord().control_volumes)}
            {
                const auto &face_mobility{pressure_field->face_mobility};
                assert(face_mobility.face_vals_axes1.rows() == first_size - 1ll);
                assert(face_mobility.face_vals_axes1.cols() == second_size - Grid2D_t::l_margin);
                assert(face_mobility.face_vals_axes2.rows() == first_size);
                assert(face_mobility.face_vals_axes2.cols() == second_size - Grid2D_t::l_margin - 1ll);

                for (auto col{Grid2D_t::l_margin + 1ll}, count{0ll}; col < second_size; ++col, ++count)
                {
                    for (auto row{0ll}; row < mobility_factor.rows(); ++row)
                        assert(!std::isnan(mobility_factor(row, count)) && !std::isinf(mobility_factor(row, count)));
                }
            }

            /// @brief set the flow field at the next time moment
            /// @param t current time moment
            /// @param t_step step to the next time moment
            void set_flow_field(
                double t, RealType t_step)
            {
                const auto &mob_factor{mobility_factor};

                Base::set_flow_field(
                    t, t_step,
                    [&mob_factor](const ptrdiff_t count, const ptrdiff_t col, const auto &P)
                    {
                        return (mob_factor.col(count) * (P.col(col - 1ll) - P.col(col))).eval();
                    });
            }
            
        protected:
            const FaceValuesContainer mobility_factor;

        private:
            static auto set_mobility_factor(const auto &mobility, const auto &h)
            {
                const auto pi{std::numbers::pi_v<RealType>};
                return FaceValuesContainer{mobility.colwise() * (h * 2.0 * pi)};
            }
        };

    } // FaceProperties
} // GPN
