using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PathDependentOption.Mathematics;

namespace PathDependentOption.ExcelInterface
{
    public static class ExcelConverter
    {
        public static Matrix CreateMatrix(object[,] tab)
        {
            var M = new Matrix(tab.GetLength(0), tab.GetLength(1));
            for (int i = 0; i < tab.GetLength(0); i++)
            {
                for (int j = 0; j < tab.GetLength(1); j++)
                {
                    M[i, j] = (double)tab[i, j];
                }
            }
            return M;
        }

        public static double[,] ConvertMatrix(Matrix M)
        {
            var tab = new double[M.Get_Rows(), M.Get_Cols()];
            for (int i = 0; i < M.Get_Rows(); i++)
            {
                for (int j = 0; j < M.Get_Cols(); j++)
                {
                    tab[i, j] = M[i, j];
                }
            }
            return tab;
        }
    }
}
