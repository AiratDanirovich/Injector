#pragma once

#include <cassert>

#include <Injector/Grids/Defines.h>

namespace GPN
{
    namespace EqSolver
    {
        template <typename Coefs_t>
        struct EquationView
        {
            using Matrix_t = Eigen::SparseMatrix<RealType>; // Eigen::Block<Eigen::SparseMatrix<RealType>, 1, -1, false>;
            EquationView(Matrix_t &A,
                         RealType &rhs,
                         const ptrdiff_t diag_id,
                         const Coefs_t &neib_ids)
                : matrix{A}, rhs{rhs},
                  diag_id{diag_id},
                  neib_ids{neib_ids}
            {
                assert(matrix.cols() > 1ll);
                for (const auto id : neib_ids)
                {
                    assert(id >= 0ll);
                    assert(id < matrix.cols() * matrix.cols());
                }
            }

            Matrix_t &matrix;
            RealType &rhs;
            const ptrdiff_t diag_id;
            const Coefs_t &neib_ids;

            void set_type_I(const auto desc)
            {
                // set diagonal value = 1.0
                matrix.coeffRef(diag_id, diag_id) = desc.factor;
                // set non-diagonal values = 0.0
                for (const auto id : neib_ids)
                    matrix.coeffRef(diag_id, id) = 0.0;
                //  set rhs = val to satisfy: 1.0*T = val
                rhs = desc.value;
            }
            void add_rhs_type_II(const auto desc)
            {
                // add the given flux to the rhs,
                // if type_I BC was not applied
                // from the other face
                if (is_type_I() == false)
                    rhs += desc.value;
            }
            void set_type_III(const auto desc)
            {
                if (is_type_I() == false)
                {
                    rhs += desc.value * desc.factor;
                    matrix.coeffRef(diag_id, diag_id) += desc.factor;
                }
            }

        protected:
            bool is_type_I() const
            { // type I BC in the row, if main_diag == 1.0, and other cols == 0.0
                bool flag{matrix.coeffRef(diag_id, diag_id) == 1.0};
                for (const auto id : neib_ids)
                    flag = flag && (matrix.coeffRef(diag_id, id) == 0.0);
                return flag;
            }
        };
    } // EqSolver
} // GPN
