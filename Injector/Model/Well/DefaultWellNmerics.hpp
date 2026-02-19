#pragma once



namespace GPN
{
    namespace Wells
    {
        template <ptrdiff_t added_vars>
        struct DefaultWellNumerics
        {
             constexpr ptrdiff_t added_vars_count()
             {
                return added_vars;
             }
        };

    } // Wells

} // GPN