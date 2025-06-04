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

        struct InjectorTemperature
            : public StepPropertyGrid,
              // so far it is assumed that the rates are positive.
              // Injector
              private AssertNonNegative
        {
            InjectorTemperature(
                const StepPropertyGrid &temps)
                : StepPropertyGrid{temps},
                  AssertNonNegative{temps}
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
        History(const Logs::InjectorRate &rates,
                const Logs::InjectorTemperature &temps)
            : rates{rates},
              temps{temps},
              time_steps{rates.grid.dual_steps},
              time_moments(time_steps.size() + 1ll, 0.0)
        {
            std::partial_sum(
                time_steps.cbegin(),
                time_steps.cend(),
                time_moments.begin() + 1ull);
        }

        const Logs::InjectorRate rates;
        const Logs::InjectorTemperature temps;
        const DualStepsContainer time_steps;
        std::vector<double> time_moments;
    };

    struct HistoryFactory
    {
        static auto create(
            const auto &time_steps,
            const auto &rates,
            const auto &temps)
        {
            const auto time{Grids::TemporalGridDual{
                Grids::GridDualStencils{
                    DualStepsContainer{time_steps}},
                CoordinateTypes::Time{}}};
            return History{
                Logs::InjectorRate{Logs::StepPropertyGrid{
                    rates, time}},
                Logs::InjectorTemperature{Logs::StepPropertyGrid{
                    temps, time}}};
        }
    };
} // GPN
