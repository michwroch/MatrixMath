using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public class Macierz
    {
        #region Definitions
        public double[,] Matrix;
        public int rows;
        public int cols;

        public Macierz(int iRows, int iCols)
        {
            rows = iRows;
            cols = iCols;
            Matrix = new double[rows, cols];
        }
        
        public double this[int iRow, int iCol]
        {
            get { return Matrix[iRow, iCol]; }
            set { Matrix[iRow, iCol] = value; }
        }

        public int GetLength(int i)
        {
            int outp = 0;
            if (i == 0)
                outp = rows;
            else
                outp = cols;
            return outp;
        }
        #endregion

        #region Print
        public override string ToString()
        {
            string output = string.Empty;
            int max1 = Matrix.GetLength(0);
            int max2 = Matrix.GetLength(1);
            for (int i = 0; i < max1; i++)
            {
                for (int j = 0; j < max2; j++)
                {
                    if (j != max2 - 1)
                        output = output + Matrix[i, j].ToString() + " ";
                    else
                        output = output + Matrix[i, j].ToString() + "\n";
                }
            }
            return output;
        }
        #endregion

        #region SUBTRACTION
        //SUBTRACTION
        public static Macierz Roznica(Macierz A, Macierz B)
        {
            Macierz output = new Macierz(A.rows, A.cols);

            if (A.cols == B.cols && A.rows == B.rows)
            {
                for (int i = 0; i < A.rows; i++)
                {
                    for (int j = 0; j < A.cols; j++)
                    {
                        output[i, j] = A[i, j] - B[i, j];
                    }
                }
            }
            else
            {

            }
            return output;
        }
        public static Macierz operator -(Macierz m1, Macierz m2)
        { return Macierz.Roznica(m1, m2); }
        #endregion
        #region ADDITION
        //ADDITION
        public static Macierz Suma(Macierz A, Macierz B)
        {
            Macierz output = new Macierz(A.GetLength(0), A.GetLength(1));

            if (A.GetLength(0) == B.GetLength(0) && A.GetLength(1) == B.GetLength(1))
            {
                for (int i = 0; i < A.GetLength(0); i++)
                {
                    for (int j = 0; j < A.GetLength(1); j++)
                    {
                        output[i, j] = A[i, j] + B[i, j];
                    }
                }
            }
            else
            {

            }
            return output;
        }
        public static Macierz operator +(Macierz m1, Macierz m2)
        { return Macierz.Suma(m1, m2); }

        #endregion
        #region Constant
        //Constant
        public static Macierz Stala(Macierz A, double B)
        {
            Macierz output = new Macierz(A.GetLength(0), A.GetLength(1));

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    output[i, j] = A[i, j] * B;
                }
            }

            return output;
        }
        #endregion
        #region MULTIPLICATION
        //MULTIPLICATION
        public static Macierz Iloczyn(Macierz A, Macierz B)
        {
            int maxA1 = A.GetLength(0);
            int maxA2 = A.GetLength(1);
            int maxB1 = B.GetLength(0);
            int maxB2 = B.GetLength(1);

            Macierz output = new Macierz(maxA1, maxB2);
            if (maxA2 == maxB1)
            {
                for (int i = 0; i < maxA1; ++i)
                {
                    for (int j = 0; j < maxB2; ++j)
                    {
                        for (int k = 0; k < maxA2; ++k)
                        {
                            output[i, j] += A[i, k] * B[k, j];
                        }
                    }
                }
            }
            else
            {

            }

            return output;
        }
        public static Macierz operator *(Macierz m1, Macierz m2)
        { return Macierz.Iloczyn(m1, m2); }
        public static Macierz operator *(Macierz m2, double m1)
        { return Macierz.Stala(m2, m1); }
        #endregion
        #region TRANSPOSITION
        //TRANSPOSITION
        public static Macierz T(Macierz A)
        {
            Macierz output = new Macierz(A.GetLength(1), A.GetLength(0));

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    output[j, i] = A[i, j];
                }
            }

            return output;
        }
        public static Macierz operator ~(Macierz m1)
        { return Macierz.T(m1); }
        #endregion

        #region NULLMATRIX 
        //Macierz zerowa AxB
        public static Macierz Zerowa(int A, int B)
        {
            Macierz output = new Macierz(A, B);
            return output;
        }
        #endregion
        #region DIAGONALEMPTY
        //Macierz jednostkowa
        public static Macierz I(int A)
        {
            Macierz output = new Macierz(A, A);


            for (int i = 0; i < A; i++)
            {
                for (int j = 0; j < A; j++)
                {
                    if (i == j)
                    {
                        output[i, j] = 1;
                    }
                    else
                    {
                        output[i, j] = 0;
                    }
                }
            }

            return output;
        }
        #endregion

        #region CHOLESK'y
        //Dekompozycja Cholesky'ego 
        public static Macierz CholeskyD(Macierz A)
        {
            int maxA1 = A.GetLength(0);
            int maxA2 = A.GetLength(1);

            Macierz L = new Macierz(maxA1, maxA1);
            L[0, 0] = Math.Sqrt(A[0, 0]);

            if (maxA2 == maxA1)
            {
                double sk = 0.0;
                double sj = 0.0;

                for (int i = 0; i < maxA1; i++)
                {
                    //
                    sk = 0.0;
                    //
                    for (int k = 0; k < i; k++)
                    {
                        //
                        sk += Math.Pow(L[i, k], 2.0);
                    }
                    //
                    L[i, i] = Math.Sqrt(A[i, i] - sk);

                    //
                    for (int j = i + 1; j < maxA1; j++)
                    {
                        //
                        sj = 0.0;
                        //
                        for (int k = 0; k < i; k++)
                        {
                            //
                            sj += L[i, k] * L[j, k];
                        }
                        //
                        L[j, i] = (A[j, i] - sj) / L[i, i];
                    }
                    //
                }
            }
            else
            {

            }

            return L;
        }

        public static Macierz CholeskyG(Macierz A)
        {
            int maxA1 = A.GetLength(0);
            int maxA2 = A.GetLength(1);

            Macierz L = new Macierz(maxA1, maxA1);
            L = CholeskyD(A);
            Macierz LG = new Macierz(maxA1, maxA1);

            if (maxA2 == maxA1)
            {
                for (int i = 0; i < maxA1; i++)
                {
                    for (int j = 0; j < maxA1; j++)
                    {
                        if (i != j)
                        {
                            LG[i, j] = L[j, i];
                        }
                        else
                        {
                            LG[i, j] = L[i, j];
                        }
                    }
                }
            }
            else
            {

            }

            return LG;
        }
        #endregion
        #region EQ System
        //Rozwiąż układ równań
        public static Macierz Rownanie(Macierz A, Macierz B)
        {
            int max = B.GetLength(0);
            Macierz output1 = new Macierz(max, 1);
            Macierz output = new Macierz(max, 1);

            int maxA = A.GetLength(0);
            int maxA2 = A.GetLength(1);
            Macierz G = new Macierz(maxA, maxA2);
            Macierz D = new Macierz(maxA, maxA2);
            G = CholeskyG(A);
            D = CholeskyD(A);

            Macierz input = new Macierz(maxA, maxA2);
            double s = 0;

            for (int i = 0; i < maxA; i++)
            {
                s = 0;
                for (int k = 0; k < maxA; k++)
                {
                    s += D[i, k] * output1[k, 0];
                }
                output1[i, 0] = ((B[i, 0] - s) / D[i, i]);
            }

            for (int i = maxA - 1; i >= 0; i--)
            {
                s = 0;
                for (int k = maxA - 1; k >= 0; k--)
                {
                    s += G[i, k] * output[k, 0];
                }
                output[i, 0] = ((output1[i, 0] - s) / G[i, i]);
            }

            return output;
        }
        public static Macierz operator /(Macierz m1, Macierz m2)
        { return Macierz.Rownanie(m2, m1); }
        #endregion

        private static Random random;
        private static object syncObj = new object();
        private static int GenerateRandomNumber(int max)
        {
            lock (syncObj)
            {
                if (random == null)
                    random = new Random(); // Or exception...
                return random.Next(max);
            }
        }
        #region RANDOM MATRIX  
        //Macierz lodsowa
        public void Losowa()
        {
            Random Rand = new Random();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    this[i, j] = GenerateRandomNumber(1000);
                }
            }
        }
        #endregion

        #region ROUNDING CELLS

        //Zaokrąglanie 
        public static Macierz Zaok_macierz(Macierz input, int r)
        {
            int w = input.GetLength(0);
            int k = input.GetLength(1);
            Macierz output = new Macierz(w, k);
            for (int a = 0; a < w; a++)
            {
                for (int b = 0; b < k; b++)
                {
                    output[a, b] = Math.Round(input[a, b], r);
                }
            }
            return output;
        }

        #endregion

    }

    public static class MacierzHelper
    {
        #region Conversion from, to
        //Konwersja na double
        public static Macierz ToMacierz(this int[,] input)
        {
            Macierz output = new Macierz(input.GetLength(0), input.GetLength(1));
            for (int i = 0; i < input.GetLength(0); i++)
            {
                for (int j = 0; j < input.GetLength(1); j++)
                {
                    output[i, j] = input[i, j];
                }
            }
            return output;
        }

        public static double[] ToArray(this Macierz M)
        {
            double[] Out = new double[M.rows];

            for (int i = 0; i < M.rows; i++)
            {
                Out[i] = M[i, 0];
            }

            return Out;
        }

        public static Macierz ToMacierz(this double[] M)
        {
            Macierz Out = new Macierz(M.GetLength(0), 1);

            for (int i = 0; i < M.GetLength(0); i++)
            {
                Out[i, 0] = M[i];
            }

            return Out;
        }

        #endregion
    }
}
