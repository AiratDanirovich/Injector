#pragma once

#include <Injector/Declarations.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Grids
    {
        /// @brief Concrete grid along radial direction of cylinder system of coordinates
        struct RGrid : public AxesGrid<CoordinateTypes::R_CylCoord>{};
        /// @brief Concrete grid along Z direction of Cartesian system of coordinates
        struct ZGrid : public AxesGrid<CoordinateTypes::Z>{};

        /// @brief Concrete grid along X direction of Cartesian system of coordinates
        struct XGrid : public AxesGrid<CoordinateTypes::X>{};
        /// @brief Concrete grid along Y direction of Cartesian system of coordinates
        struct YGrid : public AxesGrid<CoordinateTypes::Y>{};
    } // Grids

} // GPN

