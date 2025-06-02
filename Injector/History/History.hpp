#pragma once
#include <numeric>

#include <Injector/Properties/Logs.hpp>
#include <Injector/History/TemporalGrid.hpp>

namespace GPN
{
    namespace Logs
    {
        struct InjectorRate
            : public StepPropertyGrid,
              // so far it is assumed that the rates are positive.
              // Injector
              private AssertNonNegative
        {
            InjectorRate(
                const StepPropertyGrid &rates)
                : StepPropertyGrid{rates},
                  AssertNonNegative{rates}
            {
            }
        };

        struct BottomholePressure
            : public StepPropertyGrid,
              private AssertNonNegative
        {
            BottomholePressure(
                const StepPropertyGrid &pressure)
                : StepPropertyGrid{pressure},
                  AssertNonNegative{pressure}
            {
            }
        };
    } // Logs

    struct History
    {
        History(const Logs::InjectorRate &rates)
            : rates{rates},
              time_steps{rates.grid.dual_steps},
              time_moments(time_steps.size() + 1ll, 0.0)
        {
            std::partial_sum(
                time_steps.cbegin(),
                time_steps.cend(),
                time_moments.begin() + 1ull);
        }

        const Logs::InjectorRate rates;
        const DualStepsContainer time_steps;
        std::vector<double> time_moments;
    };

    struct HistoryFactory
    {
        static auto create(
            const auto &time_steps,
            const auto &rates)
        {
            const auto temp0{DualStepsContainer{time_steps}};
            const auto temp{Grids::GridDualStencils{temp0}};
            const auto temp2{Grids::TemporalGridDual{
                        temp, CoordinateTypes::Time{}}};
            return History{
                Logs::InjectorRate{Logs::StepPropertyGrid{
                    rates, temp2
                    }}};
        }
    };
} // GPN
