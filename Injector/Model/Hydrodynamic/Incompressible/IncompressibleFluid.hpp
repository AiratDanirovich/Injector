#pragma once

#include <numbers>
#include <cassert>

#include <Injector/Solver/State2D.hpp>

#include <Injector/Model/Phases/PhaseProperties.hpp>
#include <Injector/Properties/Logs.hpp>
#include <Injector/Properties/PhysicalField.hpp>

namespace GPN
{
    namespace Hydrodynamic
    {
        template <typename Grid2D_t, typename Fluid_t, typename Well_t>
        struct IncompressibleFluidField
        {
            IncompressibleFluidField(
                const RealType start_time,
                const Fluid_t &fluid,
                const Logs::Permeability &permeability,
                const Logs::ExternalPressure &ext_pressure,
                const Well_t &well,
                const cptr<Grid2D_t> grid2D)
                : fluid{fluid},
                  permeability{permeability},
                  ext_pressure{ext_pressure},
                  P_ext{Properties::FieldFactory::create(ext_pressure, grid2D)},
                  thickness_log{grid2D->first_coord.control_volumes},
                  well{well},
                  grid2D{grid2D},
                  auxillary_term{set_auxillary_term(fluid, ext_pressure, permeability, grid2D)}
            {
            }

            void set_pressure_field(const StepPropertyContainer &RFP)
            {
                // calculate current pressure as if there is no well and cement-sandwich
                auto v{((auxillary_term.colwise() * RFP).colwise() + ext_pressure.log_vals).eval()};
                // set pressure in cement, sandwich and fluid
                // equal to the pressure in the first rocks cell
                v.leftCols(3ll).colwise() = v.col(3ll);
                P = std::make_shared<Properties::Pressure<Grid2D_t>>(
                    v,
                    grid2D);
            }

            const Properties::Pressure<Grid2D_t> &current_pressure() const
            {
                return *P;
            }

        private:
            std::shared_ptr<Properties::Pressure<Grid2D_t>> P;

        public:
            const Properties::Pressure<Grid2D_t> P_ext;
            const Fluid_t& fluid;
            const Well_t& well;
            const cptr<Grid2D_t> grid2D;
            const ControlVolumesContainer &thickness_log;
            const Logs::Permeability& permeability;
            const Logs::ExternalPressure& ext_pressure;

            const CellNodesContainer2D auxillary_term;

        private:
            static CellNodesContainer2D set_auxillary_term(
                const auto &fluid, const auto &ext_pressure, const auto &permeability, const auto grid2D)
            {
                const auto thickness_log{grid2D->first_coord.control_volumes};
                const auto &r_grid{grid2D->second_coord};
                const auto r_max{r_grid.dual_back()};
                const auto pi{std::numbers::pi_v<RealType>};

                const StepPropertyContainer temp1{(-fluid.viscosity) * (r_grid.mesh_nodes / r_max).log()};
                StepPropertyContainer temp2{1.0 / ((2.0 * pi) * thickness_log * permeability.log_vals)};

                for (auto row{0ll}; row < permeability.size(); ++row)
                    if (permeability(row) == 0.0)
                        temp2(row) = ext_pressure(row);

                const CellNodesContainer2D temp{(temp2.matrix() * temp1.transpose().matrix()).array()};

                assert(temp.rows() == grid2D->first_coord.mesh_size());
                assert(temp.cols() == grid2D->second_coord.mesh_size());

                return temp;
            }
        };
    } // Hydrodynamic
} // GPN