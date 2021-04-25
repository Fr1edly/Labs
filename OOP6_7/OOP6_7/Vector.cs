using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOP6_7
{
    class Vector : Matrix
    {
        public Vector(int r=3) : base(r, 1)
        {
            mtr_id = IDS;
            Console.WriteLine("Вызвался конструктор и создалася вектор  " + mtr_id);
        }

        public Vector(int r, double[] arr) : base(r, 1, arr)
        {
            mtr_id = IDS;
            Console.WriteLine("Вызвался конструктор и создалася вектор  " + mtr_id);
        }

        public Vector(int r, Func<int, double> genVec) : base(r, 1)
        {
            for (int i = 0; i < r; i++)
                mtr[i] = genVec(i);
            mtr_id = IDS;
            Console.WriteLine("Вызвался конструктор и создалася вектор  " + mtr_id);
        }
        public Vector(Matrix other) : base(other)
        {
            if (cols > 1)
                throw new Exception("Ошибка при копировании");
            mtr_id = IDS;
            Console.WriteLine("Вызвался конструктор и скопировался вектор " + mtr_id);
        }

        ~Vector()
        {
            Console.WriteLine("Вызвался деструктор и удалился вектор " + mtr_id);
        }

        public static Vector operator +(Vector A, Vector B)
        {
            if (A.CheckForSum(B) && !A.isNULL() && !B.isNULL())
            {
                Vector temp = new Vector(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] + B.mtr[i];
                }
                return temp;
            }
            else
                throw new ArgumentException("Ошибка! Вектора с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя складывать!");
        }

        public static Vector operator -(Vector A, Vector B)
        {
            if (A.CheckForSum(B) && !A.isNULL() && !B.isNULL())
            {
                Vector temp = new Vector(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] - B.mtr[i];
                }
                return temp;
            }
            else
                throw new ArgumentException("Ошибка! Вектора с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя вычитать!");
        }

        public static Vector operator *(Vector A, double k)
        {
            if (!A.isNULL())
            {
                Vector temp = new Vector(A);
                for (int i = 0; i < (A.rows * A.cols); i++)
                {
                    temp.mtr[i] = A.mtr[i] *= k;
                }
                return temp;
            }
            else throw new ArgumentException("Ошибка! Вектор с номером " + A.mtr_id + " пуст!");
        }
        public static Matrix operator *(Vector A, Vector B)
        {
            if (A.CheckForMult(B) && !A.isNULL() && !B.isNULL())
            {
                Matrix C = new Matrix(A.rows, B.cols);
                for (uint i = 0; i < A.rows; i++)
                {
                    for (uint j = 0; j < B.cols; j++)
                    {
                        C[i, j] = 0;
                        for (uint k = 0; k < A.cols; k++)
                            C[i, j] += A[i, k] * B[k, j];
                    }
                }
                return C;
            }
            else
                throw new ArgumentException("Ошибка! Вектора с номерами " + A.mtr_id + " и " + B.mtr_id + " нельзя умножать!");
        }

        public Vector Transponse()
        {
            int temp = cols;
            cols = rows;
            rows = temp;
            return this;
        }
        public Vector Mult_on_scal(double scal)
        {
            if (!isNULL())
            {
                for (int i = 0; i < rows; i++)               
                        mtr[i] *= scal;
                return this;
            }
            else throw new ArgumentException("Вектор " + mtr_id + " Пуст");
        }

        public Vector Nomralization()
        {
            if (!isNULL())
            {
                double module = Module();
                for (int i = 0; i < rows; i++)
                    mtr[i] /= module;
                return this;
            }
            else throw new ArgumentException("Вектор " + mtr_id + " Пуст");
        }

        public double Module()
        {
            if (!isNULL())
            {
                double result = 0;
                for (int i = 0; i < rows; i++)
                    result += mtr[i] * mtr[i];

                return Math.Sqrt(result);
            }
            else throw new ArgumentException("Вектор " + mtr_id + " Пуст");
        }

        public Vector VecMult(Vector A, Vector B)
        {
            if (A.rows > 3 || B.rows > 3) throw new Exception("Вектор больше нужного размера");
            else
            {
                Vector C = new Vector();
                C.mtr[0] = (A.mtr[1] * B.mtr[2]) - (A.mtr[2] * B.mtr[1]);
                C.mtr[1] = (A.mtr[2] * B.mtr[0]) - (A.mtr[0] * B.mtr[2]);
                C.mtr[2] = (A.mtr[0] * B.mtr[1]) - (A.mtr[1] * B.mtr[0]);

                C.cols = 1;
                C.rows = 3;
                return C;
            }
        }

        public double Cos(Vector A, Vector B)
        {
            return MultScal(A, B) / (A.Module() * B.Module());
        }

        public double Sin(Vector A, Vector B)
        {
            return 1 - Math.Pow(Cos(A, B), 2);
        }

        public double MultScal(Vector A, Vector B)
        {
            if (!A.CheckForSum(B) && !A.isNULL() && !B.isNULL())
            {
                double res = 0;
                for (int i = 0; i < A.rows; i++)           
                        res = A.mtr[i] * B.mtr[i];
                return res;
            }
            else throw new Exception("Ошибка!");
        }

        public double Angle(Vector A, Vector B)
        {
            return Math.Acos(MultScal(A, B) / (A.Module() * B.Module())) * 180 / Math.PI;
        }
    }
}
