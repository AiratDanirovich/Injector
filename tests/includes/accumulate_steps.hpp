#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

#include <Injector/Grids/Defines.h>
#include "transfer_to_vector.hpp"

using VR = std::vector<GPN::RealType>;

VR accumulate_steps(const VR &depth_intervals)
{
  VR depth_stencils(depth_intervals.size() + 1ull, 0.0);
  std::partial_sum(
      depth_intervals.cbegin(),
      depth_intervals.cend(),
      depth_stencils.begin() + 1ull);

  assert(depth_stencils[0ull] == 0.0);
  for (auto i{1ull}; i < depth_stencils.size(); ++i)
    assert(depth_stencils[i] == depth_stencils[i - 1] + depth_intervals[i - 1ll]);

  return depth_stencils;
}

VR accumulate_steps(const Eigen::ArrayX<GPN::RealType> &depth_intervals)
{
  return accumulate_steps(transfer_to_vector(depth_intervals));
}