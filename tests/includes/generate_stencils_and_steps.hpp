#pragma once

#include <vector>

#include <Injector/Grids/Defines.h>

using VR = std::vector<GPN::RealType>;

VR generate_stencils(const VR &steps, const GPN::RealType top = 0.0)
{
    VR out;
    out.reserve(steps.size() + 1ull);
    out.push_back(top);

    for (auto i{0ull}; i < steps.size(); ++i)
        out.push_back(out.back() + steps[i]);
    return out;
}

VR generate_stencils(GPN::RealType t0, GPN::RealType t1, GPN::RealType t_step_major)
{
  auto segm_count{static_cast<size_t>(std::ceil(t1 - t0) / t_step_major)};
  double step = (t1 - t0) / (segm_count);
  VR out(segm_count + 1ll);

  for (auto i{0ull}; i < out.size(); ++i)
    out[i] = t0 + i * step;
  return out;
}

VR generate_steps(const VR &dual_nodes)
{
  VR out(dual_nodes.size() - 1ll);

  for (auto i{0ull}; i < out.size(); ++i)
    out[i] = dual_nodes[i + 1] - dual_nodes[i];
  return out;
}