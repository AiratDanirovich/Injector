#pragma once

#include <Injector/Properties/Logs.hpp>

namespace GPN
{
    struct History : public Logs::StepPropertyGrid
    {
        using StepPropertyGrid::StepPropertyGrid;
    };
} // GPN

// using DatesContainer = custom_vector<RealType>;
// using RatesContainer = custom_vector<RealType>;
// using TimeStepsContainer = custom_vector<RealType>;

// struct History
//     {
//         History(
//             const DatesContainer &time_step_stencils_,
//             const RatesContainer& rates,
//             const DatesContainer& time_moments_) noexcept
//             : time_step_stencils{time_step_stencils_}
//             , time_stencils{make_time_stencils(time_step_stencils_)}
//             , rates{rates}
//             , time_moments{time_moments}
//             , time_steps(time_moments_.size()-1)
//         {
//             for(size_t id{0ull}; id < time_steps.size(); ++id)
//                 time_steps(id) = time_moments(id+1) - time_moments(id);
//         }

//         template <typename RefinementPolicy>
//         History(
//             const DatesContainer &time_step_stencils,
//             const RatesContainer& rates,
//             const RefinementPolicy& refiner) noexcept
//             : History{time_stencils, rates, refiner.refine(time_moments)}
//         {}

//     protected:
//         // time moments when the rate is changed
//         DatesContainer time_step_stencils;
//         DatesContainer time_stencils;
//         RatesContainer rates;

//         DatesContainer time_moments;
//         TimeStepsContainer time_steps;

//     protected:
//         static auto make_time_stencils(const DatesContainer &time_step_stencils)
//         {
//             DatesContainer out(time_step_stencils.size()+1);
//             out(0) = 0;
//             for(size_t id{1ull}; id < time_step_stencils.size(); ++id)
//                 out(id) = out(id-1) + time_step_stencils(id-1);
//             return out;
//         }
//     };