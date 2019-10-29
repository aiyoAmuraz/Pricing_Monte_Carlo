using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExcelDna.Integration;
using PathDependentOption.Mathematics;

namespace PathDependentOption.ExcelInterface
{
    public static class ExcelFunction
    {
        [ExcelFunction(Description = "Make sure the matrix is symmetric")]
        public static bool IsSymmetric(object[,] tab)
        {
            var M = ExcelConverter.CreateMatrix(tab);
            return M.IsSymmetric();
        }

        [ExcelFunction(Description = "Get the lower triangular matrix from the Cholesky decomposition")]
        public static double[,] CholeskyDecomposition(object[,] tab)
        {
            var M = ExcelConverter.CreateMatrix(tab);
            return ExcelConverter.ConvertMatrix(M.CholeskyDecomposition());
        }
    }
}
