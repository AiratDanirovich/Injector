#pragma once
#include <numeric>
#include <limits>
#include <cassert>

#include <Injector/History/InjectorRegimes.hpp>

#include <Injector/Properties/Logs.hpp>
#include <Injector/History/TemporalGrid.hpp>
#include <Injector/History/InjectorRegimes.hpp>

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
        struct TimeStep : public SomeProperty
        {
        };
        struct StartTime : public SomeProperty
        {
        };

        struct Record
        {
            Record(
                const Pressure pressure,
                const Rate rate,
                const InjectorRegimes::Type type,
                const StartTime start_time,
                const TimeStep time_step)
                : pressure{pressure}, rate{rate}, type{type},
                  start_time{start_time}, time_step{time_step},
                  end_time{start_time+time_step},
                  mid_time{start_time + time_step/2.0}
            {
            }
            const RealType rate;
            const RealType pressure;
            const InjectorRegimes::Type type;
            const RealType start_time, time_step, end_time, mid_time;
        };

        History(const Logs::InjectorRate &rates,
                const Logs::SurfacePressure &pressures,
                const Logs::InjectorTemperature &temps,
                const std::vector<InjectorRegimes::Type> &regimes)
            : rates{rates},
              pressures{pressures},
              temps{temps},
              time_steps{rates.grid.dual_steps},
              time_moments{set_time_moments(rates.grid.dual_steps)},
              regimes{regimes},
              pos{-1ll}
        {
            assert(rates.size() == time_steps.size());
            assert(pressures.size() == time_steps.size());
            assert(temps.size() == time_steps.size());
            assert(regimes.size() == time_steps.size());

            for (auto i{0ll}; i < rates.size(); ++i)
            {
                if (regimes[i] == InjectorRegimes::FixedBottomHolePressure)
                {
                    assert(std::isnan(rates.log_vals(i)));
                    assert(pressures.log_vals(i) > 0.0);
                }
                else if (regimes[i] == InjectorRegimes::FixedRate)
                {
                    assert(std::isnan(pressures.log_vals(i)));
                    //    assert(rates.log_vals(i) >= 0.0);
                }
                else
                    assert(false && "Wrong injector regime!");
            }
        }

        void reset()
        {
            pos = -1ll;
        }

        void advance()
        {
            ++pos;
        }
        const auto size() const
        {
            return rates.size();
        }
        const auto record_id() const
        {
            return pos;
        }

        const auto temperature() const
        {
            assert(pos >= 0ll);
            assert(pos < (ptrdiff_t)size());
            return temps(pos);
        }
        const auto pressure() const
        {
            assert(pos >= 0ll);
            assert(pos < (ptrdiff_t)size());
            return pressures(pos);
        }

        const auto regime() const
        {
            assert(pos >= 0ll);
            assert(pos < (ptrdiff_t)size());
            return regimes[pos];
        }
        const auto rate() const
        {
            assert(pos >= 0ll);
            assert(pos < (ptrdiff_t)size());
            return rates(pos);
        }

        const auto get_record(auto idx) const
        {
            assert(idx >= 0ll);
            assert(idx < (ptrdiff_t)size());
            return Record{
                Pressure{pressures(idx)}, 
                Rate{rates(idx)}, 
                regimes[idx],
                StartTime{time_moments[idx]},
                TimeStep{time_steps(idx)}};
        }

        const auto get_current_record() const
        {
            assert(pos >= 0ll);
            assert(pos < (ptrdiff_t)size());
            return get_record(pos);
        }

        const Logs::InjectorRate rates;
        const Logs::SurfacePressure pressures;
        const std::vector<InjectorRegimes::Type> regimes;
        const Logs::InjectorTemperature temps;
        const DualStepsContainer time_steps;
        const std::vector<double> time_moments;

    private:
        ptrdiff_t pos;

        std::vector<double> set_time_moments(const DualStepsContainer &time_steps)
        {
            std::vector<double> time_moments(time_steps.size() + 1ll, 0.0);
            std::partial_sum(
                time_steps.cbegin(),
                time_steps.cend(),
                time_moments.begin() + 1ull);
            return time_moments;
        }
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
            const std::vector<RealType> p(rates.size(),
                                          std::numeric_limits<double>::quiet_NaN());

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
                InjectorRegimes::FixedBottomHolePressure);

            // nan-valued rate.
            // can be calculated during the simulation
            const std::vector<RealType> q(pressure.size(),
                                          std::numeric_limits<double>::quiet_NaN());

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
