#include <iostream>
#include <math.h>
#include <iomanip>

#include "slave.h"

namespace slave {

	double* orto(double** slay, double* B_mass, const size_t len) {

		double* x = new double[len](); //вектор для ответов
		double **U = new double*[len]; // буферный массив для Ui
		for (int i = 0; i < len; i++) U[i] = new double[len]();

		double** Z = new double* [len]; // решение СЛАУ
		for (int i = 0; i < len; i++) Z[i] = new double[len]();

		double *L = new double[len](); //L с N=3 значениями  - это массив для  ..хранения подкоренного скалярного произведения(ui, ui), т.е. длина вектора
		double* TMP = new double[len];
		double* H = new double[len];   //для B_mass значений

		//Начальная иттерация:
		for (int i = 0; i < len; i++) {
			U[0][i] = slay[0][i]; 
		}

		L[0] = sqrt(scalar(U[0], U[0], len));

		for (int i = 0; i < len; i++) {
			Z[0][i] = U[0][i] * (1 / L[0]);
		}
		double* T = new double[len];
		T[0] = B_mass[0] / L[0];


///последующие иттерации
		int count = 1; ///Счетчик переменных для которых уже найдены первые значения 
		int n = 0;
		double t;
		
		do {
			double temp_h; 
			temp_h = 0;

			for (int i = 0; i < len; i++) TMP[i] = 0;
			///U_i
			for (int j = 0; j <= count - 1; j++) {
				t = scalar(slay[count], Z[j], len);	///Получили c[i,j] 
				for (int i = 0; i < len; i++) {
					TMP[i] = TMP[i] + (t * Z[j][i]); ///Получили c[i,j]*z[k] // (a[i],z[k])*z[k]
				}
				temp_h += t * T[j]; /// (a[i],z[k])*(-b(f)[i]/len_vect)
			}
			///Нашли U
			for (int i = 0; i < len; i++)
				U[count][i] = slay[count][i] - TMP[i];///  a[i] - (a[i],z[k])*z[k]
			///Нашли L len_vect
			L[count] = sqrt(scalar(U[count], U[count], len));
			///Нашли H
			T[count] = (B_mass[count] - temp_h) / (L[count]);
			///Z_[i]
			for (int i = 0; i < len; i++) {
				Z[count][i] = (1 / L[count] * (U[count][i]));
			}
			n++;
			count++;
		} while (count < len);
		//
		for (int i = 0; i < len; i++) 
			for (int j = 0; j < len; j++)
				x[i] += Z[j][i] * T[j];
	
#pragma region PURGE
		for (int i = 0; i < len; i++) delete[] U[i];
		delete[] U;
		for (int i = 0; i < len; i++) delete[] Z[i];
		delete[] Z;
		delete[] L;
		delete[] TMP;
		delete[] H;
#pragma endregion
		
		return x;
	}
	/// скалярное произведение
	double scalar(double* fobj, double* sobj, double len) {
		double ans = 0;
		for (int j = 0; j < len; j++) {
			ans += fobj[j] * sobj[j];
		}
		return ans;
	}

	double* iter(double** slay, double* B_mass, const size_t len, double eps) {
		double* x, * x0, * E, max;
		x = new double[len]();
		x0 = new double[len]();
		E = new double[len]();
		size_t count = 0;
		do {
			for (size_t i = 0; i < len; i++) {
				for (int j = 0; j < len; j++)
					x[i] += slay[i][j] * x0[j];
				x[i] += B_mass[i];
				E[i] = fabs(x[i] - x0[i]);
			}
			max = 0;
			int i;
			for (i = 0; i < len; i++) {
				if (max < E[i]) max = E[i];
				x0[i] = x[i];
			}
			count++;
		} while (max > eps);
		std::cout << "Количество итераций: " << count << std::endl << std::endl;
		return x;
	}

	double* miter(double** a, double* y, const size_t n)
	{
		double* res = new double[n];
		int i, j;


		for (i = 0; i < n; i++)
		{
			res[i] = y[i] / a[i][i];
		}

		double eps = 0.0001;
		double* Xn = new double[n];

		do {
			for (i = 0; i < n; i++) {
				Xn[i] = y[i] / a[i][i];
				for (j = 0; j < n; j++) {
					if (i == j)
						continue;
					else {
						Xn[i] -= a[i][j] / a[i][i] * res[j];
					}
				}
			}

			bool flag = true;
			for (i = 0; i < n - 1; i++) {
				if (abs(Xn[i] - res[i]) > eps) {
					flag = false;
					break;
				}
			}

			for (i = 0; i < n; i++) {
				res[i] = Xn[i];
			}

			if (flag)
				break;
		} while (1);

		return res;
	}
}