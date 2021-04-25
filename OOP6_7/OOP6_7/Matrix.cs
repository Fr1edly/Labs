using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOP6_7
{
    class Matrix
    {

        protected double[] mtr;
        protected int rows, cols;
        protected uint mtr_id;
        protected static uint IDS = 0;


        public Matrix(int r, int c, Func<uint, uint, double> genMtr)
        {
            if (r == 0 || c == 0)
            {
                mtr_id = ++IDS;
                this.cols = 0;
                this.rows = 0;
                mtr = null;
                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
            else
            {
                mtr_id = ++IDS;
                this.cols = c;
                this.rows = r;
                mtr = new double[r * c];
                for (uint i = 0; i < r; i++)
                {
                    for (uint j = 0; j < c; j++)
                        mtr[i * c + j] = genMtr(i, j);
                }
                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
        }

        public Matrix(int r, int c)
        {
            if (r == 0)
            {
                mtr_id = ++IDS;
                this.rows = 0;
                this.cols = 0;
                mtr = null;

                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
            else
            {
                mtr_id = ++IDS;
                this.rows = r;
                this.cols = c;
                mtr = new double[r * c];
                Random rand = new Random();
                for (int i = 0; i < r * c; i++)
                    mtr[i] = rand.NextDouble();

                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
        }
        public Matrix(int s = 3)
        {
            if (s == 0)
            {
                mtr_id = ++IDS;
                rows = 0;
                cols = 0;
                mtr = null;

                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
            else
            {
                mtr_id = ++IDS;
                rows = s;
                cols = s;
                mtr = new double[s * s];
                Random rand = new Random();
                for (int i = 0; i < s * s; i++)
                    mtr[i] = rand.NextDouble();
                Console.WriteLine("Вызвался конструктор и создалась матрица " + mtr_id);
            }
        }
        public Matrix(int r, int c, double[] arr)
        {
            if (r != 0 && c != 0)
            {
                this.rows = r;
                this.cols = c;
                mtr = new double[r * c];
                for (int i = 0; i < r * c; i++)
                {
                    mtr[i] = arr[i];
                }
            }
            else Console.WriteLine("Ошибка размера матрицы");
        }

        public Matrix(Matrix other)
        {
            mtr_id = ++IDS;
            this.rows = other.rows;
            this.cols = other.cols;
            this.mtr = new double[this.rows * this.cols];

            for (int i = 0; i < (this.rows * this.cols); i++)
            {
                this.mtr[i] = other.mtr[i];
            }

            Console.WriteLine("Вызвался конструктор и скопировалась матрица " + mtr_id);
        }
        ~Matrix()
        {
            Console.WriteLine("Вызвался деструктор и удалилась матрица " + mtr_id);
            mtr_id = --IDS;
        }

        public override string ToString()
        {
            string line = "";
            if (rows == 0 || cols == 0)
            {
                line = " ";
            }
            else
            {
                for (int i = 0; i < this.rows; i++)
                {
                    for (int j = 0; j < this.cols; j++)
                    {
                        line += string.Format("{0:F2} ", this.mtr[i * this.cols + j]) + "\t";
                    }
                    line += "\n";
                }
            }
            return line;
        }

        public bool CheckForSum(Matrix a)
        {
            return ((this.rows == a.rows) && (this.cols == a.cols));
        }

        public bool CheckForMult(Matrix a)
        {
            return this.rows == a.cols;

        }
        public double MtrMax()
        {
            if (!this.isNULL())
            {
                double a = mtr[0];
                for (int i = 1; i < this.rows * this.cols; i++)
                    if (this.mtr[i] > a)
                        a = this.mtr[i];
                return a;
            }
            else
                throw new Exception("Матрица " + mtr_id + " пуста");

        }

        public double MtrMin()
        {
            if (!this.isNULL())
            {
                double a = mtr[0];
                for (int i = 1; i < this.rows * this.cols; i++)
                    if (this.mtr[i] < a)
                        a = this.mtr[i];
                return a;
            }
            else
                throw new Exception("Матрица " + mtr_id + " пуста");
        }

        public bool isNULL()
        {
            return this.mtr == null;
        }

        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A.CheckForSum(B) && !A.isNULL() && !B.isNULL())
            {
                Matrix temp = new Matrix(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] + B.mtr[i];
                }
                return temp;
            }
            else
                throw new ArgumentException("Ошибка! Матрицы с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя складывать!");
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.CheckForSum(B) && !A.isNULL() && !B.isNULL())
            {
                Matrix temp = new Matrix(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] - B.mtr[i];
                }
                return temp;
            }
            else
                throw new ArgumentException("Ошибка! Матрицы с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя вычитать!");
        }

        public static Matrix operator *(Matrix A, double k)
        {
            if (!A.isNULL())
            {
                Matrix temp = new Matrix(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] *= k;
                }
                return temp;
            }
            else throw new ArgumentException("Ошибка! Матрица с номером " + A.mtr_id + " пуста!");
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.CheckForMult(B) && !A.isNULL() && !B.isNULL())
            {
                Matrix C = new Matrix(A.rows, B.cols);
                for (int i = 0; i < A.rows * B.cols; i++)
                    C.mtr[i] = 0;
                for (int i = 0; i < A.rows; ++i)
                {
                    for (int j = 0; j < B.cols; ++j)
                    {
                        for (int k = 0; k < A.cols; ++k)
                            C.mtr[i * B.cols + j] += (A.mtr[i * A.cols + k] * B.mtr[k * B.cols + j]);
                    }
                }
                return C;
            }
            else
                throw new ArgumentException("Ошибка! Матрицы с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя умножать!");
        }

        public double this[uint r, uint c]
        {
            get
            {
                if (r < rows && c < cols)
                    return this.mtr[r * this.cols + c];
                else
                    throw new ArgumentException("Выход за границы объекта " + mtr_id);
            }
            set
            {
                if (r < rows && c < cols)
                    this.mtr[r * this.cols + c] = value;
                else
                    throw new ArgumentException("Выход за границы объекта " + mtr_id);
            }
        }
    }
}
