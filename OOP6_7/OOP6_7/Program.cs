using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OOP6_7
{
    class Program
    {
        public static double genMtr(uint i, uint j)
        {
            return i + j;
        }
        public static double genVec(int t)
        {
            return t + t;
        }
        static void Main(string[] args)
        {
            /*Matrix A = new Matrix(3,2);
            Matrix B = new Matrix(2,3);
            Matrix C= A*B;
            Console.WriteLine(A.ToString());
            Console.WriteLine(B.ToString());
            Console.WriteLine(C.ToString());
            B[1,2] = 0.2;
            Console.WriteLine(B.ToString());*/


            Vector A = new Vector();
            Vector B = new Vector();
            Console.WriteLine(A.ToString());
            Console.WriteLine(B.ToString());
            B.Transponse();
            Matrix C = A * B;
            Console.WriteLine(C.ToString());
            Console.ReadLine();
        }
    }
}
