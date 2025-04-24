#pragma once

#include <Injector/Declarations.h>
#include <Injector/Grids/Grids1D.hpp>

namespace GPN
{
    namespace Grids
    {
        struct RGrid : public AxesGrid<CoordinateTypes::R_CylCoord>{};
        struct ZGrid : public AxesGrid<CoordinateTypes::Z>{};

        struct XGrid : public AxesGrid<CoordinateTypes::X>{};
        struct YGrid : public AxesGrid<CoordinateTypes::Y>{};
    } // Grids

} // GPN

