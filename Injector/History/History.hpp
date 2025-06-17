#pragma once
#include <numeric>
#include <limits>
#include <cassert>

#include <Injector/Properties/Logs.hpp>
#include <Injector/History/TemporalGrid.hpp>

namespace GPN
{
    namespace Logs
    {
        /// @brief The rate of fluid injection
        struct InjectorRate
            : public StepPropertyGrid //,
                                      // so far it is assumed that the rates are positive.
                                      // Injector
        //  private AssertNonNegative
        {
            InjectorRate(
                const StepPropertyGrid &rates)
                : StepPropertyGrid{rates} //,
            //  AssertNonNegative{rates}
            {
            }
        };

        /// @brief Temperature of injected fluid
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

        /// @brief Surface at the to of well, P_{top}
        struct SurfacePressure
            : public StepPropertyGrid //,
        //  private AssertNonNegative
        {
            SurfacePressure(
                const StepPropertyGrid &pressure)
                : StepPropertyGrid{pressure}
            //,
            //  AssertNonNegative{pressure}
            {
            }
        };
    } // Logs

    struct InjectorRegimes
    {
        enum Type
        {
            FixedRate,
            FixedPressure
        };
    };

    struct History
    {
        struct SomeProperty
        {
            operator RealType() const { return value; }
            RealType value;
        };
        struct Pressure : public SomeProperty
        {
        };
        struct Rate : public SomeProperty
        {
        };

        struct Record
        {
            Record(Pressure pressure, Rate rate, const InjectorRegimes::Type type)
                : pressure{pressure}, rate{rate}, type{type}
            {
            }
            const RealType rate;
            const RealType pressure;
            const InjectorRegimes::Type type;
        };

        History(const Logs::InjectorRate &rates,
                const Logs::SurfacePressure &pressure,
                const Logs::InjectorTemperature &temps,
                const std::vector<InjectorRegimes::Type> &regimes)
            : rates{rates},
              pressure{pressure},
              temps{temps},
              time_steps{rates.grid.dual_steps},
              time_moments(time_steps.size() + 1ll, 0.0),
              regimes{regimes}
        {
            assert(rates.size() == time_steps.size());
            assert(pressure.size() == time_steps.size());
            assert(temps.size() == time_steps.size());
            assert(regimes.size() == time_steps.size());

            for (auto i{0ll}; i < rates.size(); ++i)
            {
                if (regimes[i] == InjectorRegimes::FixedPressure)
                {
                    assert(std::isnan(rates.log_vals(i)));
                    assert(pressure.log_vals(i) > 0.0);
                }
                else if (regimes[i] == InjectorRegimes::FixedRate)
                {
                    assert(std::isnan(pressure.log_vals(i)));
                    assert(rates.log_vals(i) > 0.0);
                }
                else
                    assert(false && "Wrong injector regime!");
            }

            std::partial_sum(
                time_steps.cbegin(),
                time_steps.cend(),
                time_moments.begin() + 1ull);
        }

        const auto get_record(auto idx) const
        {
            return Record{Pressure{pressure(idx)}, Rate{rates(idx)}, regimes[idx]};
        }

        const Logs::InjectorRate rates;
        const Logs::SurfacePressure pressure;
        const Logs::InjectorTemperature temps;
        const DualStepsContainer time_steps;
        std::vector<double> time_moments;

        const std::vector<InjectorRegimes::Type> regimes;
    };

    struct HistoryFactory
    {
        static auto createFixedRate(
            const auto &time_steps,
            const auto &rates,
            const auto &temps)
        {
            const std::vector<InjectorRegimes::Type> regimes(
                rates.size(),
                InjectorRegimes::FixedRate);

            // nan-valued pressure.
            // can be calculated during the simulation
            const std::vector<RealType> p(rates.size(), std::numeric_limits<double>::quiet_NaN());

            const auto time{Grids::TemporalGridDual{
                Grids::GridDualStencils{
                    DualStepsContainer{time_steps}},
                CoordinateTypes::Time{}}};
            return History{
                Logs::InjectorRate{Logs::StepPropertyGrid{
                    rates, time}},
                Logs::SurfacePressure{Logs::StepPropertyGrid{
                    p, time}},
                Logs::InjectorTemperature{Logs::StepPropertyGrid{
                    temps, time}},
                regimes};
        }

        static auto createFixedPressure(
            const auto &time_steps,
            const auto &pressure,
            const auto &temps)
        {
            const std::vector<InjectorRegimes::Type> regimes(
                pressure.size(),
                InjectorRegimes::FixedPressure);

            // nan-valued rate.
            // can be calculated during the simulation
            const std::vector<RealType> q(pressure.size(), std::numeric_limits<double>::quiet_NaN());

            const auto time{Grids::TemporalGridDual{
                Grids::GridDualStencils{
                    DualStepsContainer{time_steps}},
                CoordinateTypes::Time{}}};
            return History{
                Logs::InjectorRate{Logs::StepPropertyGrid{
                    q, time}},
                Logs::SurfacePressure{Logs::StepPropertyGrid{
                    pressure, time}},
                Logs::InjectorTemperature{Logs::StepPropertyGrid{
                    temps, time}},
                regimes};
        }
    };
} // GPN
