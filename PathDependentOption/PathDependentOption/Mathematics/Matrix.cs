using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PathDependentOption.Mathematics
{
    public class Matrix
    {
        int m_rows;
        int m_cols;
        private double[,] m_matrix;

        public Matrix()
        {
            m_rows = 0;
            m_cols = 0;
            m_matrix = new double[m_rows, m_cols];
        }
        public Matrix(int rows, int cols)
        {
            m_rows = rows;
            m_cols = cols;
            m_matrix = new double[m_rows, m_cols];
        }
        public double this[int rows, int cols]
        {
            get
            {
                return m_matrix[rows, cols];
            }
            set
            {
                m_matrix[rows, cols] = value;
            }
        }

        public void Display()
        {
            string val = null;
            for (int i = 0; i < m_rows; ++i)
            {
                for (int j = 0; j < m_cols; ++j)
                {
                    val += Math.Round(m_matrix[i, j], 3) + "\t";
                }
                Console.WriteLine(val);
                val = null;
            }
            Console.WriteLine();
        }
        public int Get_Rows()
        {
            return m_rows;
        }
        public int Get_Cols()
        {
            return m_cols;
        }
        public Matrix Transpose()
        {
            Matrix M = new Matrix(m_cols, m_rows);

            for (int i = 0; i < m_rows; ++i)
            {
                for (int j = 0; j < m_cols; ++j)
                {
                    M[j, i] = m_matrix[i, j];
                }
            }

            return M;
        }
        public Matrix Diagonal()
        {
            if (m_rows != m_cols)
            {
                throw new ArgumentException("Matrix must be nxn");
            }

            Matrix M = new Matrix(m_rows, m_cols);

            for (int i = 0; i < m_rows; ++i)
            {
                M[i, i] = m_matrix[i, i];
            }

            return M;
        }
        public Matrix ConvertToMatrix(double[,] array)
        {
            int rows = array.GetLength(0);
            int cols = array.GetLength(1);
            Matrix M = new Matrix(rows, cols);
            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    M[i, j] = array[i, j];
                }
            }
            return M;
        }
        public Matrix ConvertToMatrix(double[] array)
        {
            int rows = array.GetLength(0);
            Matrix M = new Matrix(rows, 1);
            for (int i = 0; i < rows; ++i)
            {
                M[i, 0] = array[i];
            }
            return M;
        }
        public double[] ConvertToVector(Matrix M)
        {
            double[] vec = new double[M.Get_Rows()];
            for (int i = 0; i < M.Get_Rows(); ++i)
            {
                vec[i] = M[i, 0];
            }
            return vec;
        }
        public bool IsSymmetric()
        {
            if (!this.IsSquare()) return false;
            for (int i = 0; i < this.m_rows; i++)
            {
                for (int j = i + 1; j < this.m_cols; j++)
                {
                    if (m_matrix[i, j] != m_matrix[j, i]) return false;
                }
            }
            return true;
        }
        public bool IsSquare()
        {
            return this.m_cols == this.m_rows ? true : false;
        }
        public Matrix CholeskyDecomposition()
        {
            if (! this.IsSymmetric()) throw new ArgumentException("Not a symmetric matrix");

            var L = new Matrix(m_rows, m_cols);
            
            for (int i = 0; i < m_rows; i++)
            {
                if (i == 0) L[i, i] = m_matrix[i, i];
                int count_rows = 0;
                for (int j = 0; j <= i; j++)
                {
                    if (i != 0)
                    {
                        if (j == 1) L[i, j] = m_matrix[1, i] / L[1, 1];
                        else if (j == i) L[i, j] = Math.Sqrt(m_matrix[i, j] - this.SumMultOfTwoLines(i,i,0,i));
                        else L[i, j] = (m_matrix[j, i] - this.SumMultOfTwoLines(i, j, 0, i)) / L[j,j];
                    }
                }
            }
            return L;
        }
        public double SumMultOfTwoLines(int row1, int row2, int n_start, int n_end)
        {
            if (n_start > n_end) throw new ArgumentException("Issue, n_start bigger than n_end");
            double res = 0;
            for (int i = n_start; i < n_end; i++)
            {
                res += m_matrix[row1, i] * m_matrix[row2, i];
            }
            return res;
        }
                          


        public static Matrix operator +(Matrix lhs, Matrix rhs)
        {
            if (lhs.m_rows != rhs.m_rows || lhs.m_cols != rhs.m_cols)
            {
                throw new ArgumentException("Index out of range");
            }

            Matrix M = new Matrix(lhs.m_rows, lhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < lhs.m_cols; ++j)
                {
                    M[i, j] = lhs[i, j] + rhs[i, j];
                }
            }

            return M;
        }
        public static Matrix operator -(Matrix lhs, Matrix rhs)
        {
            if (lhs.m_rows != rhs.m_rows || lhs.m_cols != rhs.m_cols)
            {
                throw new ArgumentException("Index out of range");
            }

            Matrix M = new Matrix(lhs.m_rows, lhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < lhs.m_cols; ++j)
                {
                    M[i, j] = lhs[i, j] - rhs[i, j];
                }
            }

            return M;
        }
        public static Matrix operator *(Matrix lhs, Matrix rhs)
        {
            if (lhs.m_cols != rhs.m_rows)
            {
                throw new ArgumentException("Number of columns of M1 must be equal to the number of rows of M2");
            }

            Matrix M = new Matrix(lhs.m_rows, rhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < rhs.m_cols; ++j)
                {
                    for (int k = 0; k < rhs.m_rows; ++k)
                    {

                        M[i, j] += lhs.m_matrix[i, k] * rhs[k, j];
                    }
                }
            }

            return M;
        }

        public static Matrix operator +(Matrix lhs, double rhs)
        {
            Matrix M = new Matrix(lhs.m_rows, lhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < lhs.m_cols; ++j)
                {
                    M[i, j] = lhs[i, j] + rhs;
                }
            }

            return M;
        }
        public static Matrix operator -(Matrix lhs, double rhs)
        {
            Matrix M = new Matrix(lhs.m_rows, lhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < lhs.m_cols; ++j)
                {
                    M[i, j] = lhs[i, j] - rhs;
                }
            }

            return M;
        }
        public static Matrix operator *(Matrix lhs, double rhs)
        {
            Matrix M = new Matrix(lhs.m_rows, lhs.m_cols);
            for (int i = 0; i < lhs.m_rows; ++i)
            {
                for (int j = 0; j < lhs.m_cols; ++j)
                {
                    M[i, j] = lhs[i, j] * rhs;
                }
            }

            return M;
        }
        public static Matrix operator +(double lhs, Matrix rhs)
        {
            Matrix M = new Matrix(rhs.m_rows, rhs.m_cols);
            for (int i = 0; i < rhs.m_rows; ++i)
            {
                for (int j = 0; j < rhs.m_cols; ++j)
                {
                    M[i, j] = rhs[i, j] + lhs;
                }
            }

            return M;
        }
        public static Matrix operator -(double lhs, Matrix rhs)
        {
            Matrix M = new Matrix(rhs.m_rows, rhs.m_cols);
            for (int i = 0; i < rhs.m_rows; ++i)
            {
                for (int j = 0; j < rhs.m_cols; ++j)
                {
                    M[i, j] = rhs[i, j] - lhs;
                }
            }

            return M;
        }
        public static Matrix operator *(double lhs, Matrix rhs)
        {
            Matrix M = new Matrix(rhs.m_rows, rhs.m_cols);
            for (int i = 0; i < rhs.m_rows; ++i)
            {
                for (int j = 0; j < rhs.m_cols; ++j)
                {
                    M[i, j] = rhs[i, j] * lhs;
                }
            }

            return M;
        }
    }
}
