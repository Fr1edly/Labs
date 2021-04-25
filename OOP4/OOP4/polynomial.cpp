#include <iostream>
#include <iomanip>
#include "polynomial.h"
#include <string>
#include <cmath>


namespace mathPolynom
{
	unsigned polynomial::IDS = 0;

	polynomial::polynomial()
	{
		pol_id = ++IDS;
		this->size = 0;
		this->coef = new double[size + 1];
		this->coef[0] = 0.0;
#ifdef _DEBUG
		std::cout << "The polynomial " << pol_id << " created" << std::endl;
#endif
	}

	polynomial::polynomial(int S)
	{
		pol_id = ++IDS;
		this->size = S;
		this->coef = new double[size]();
#ifdef _DEBUG
		std::cout << "The polynomial " << pol_id << " created" << std::endl;
#endif
	}

	polynomial::polynomial(int S, double* C)
	{
		if (S == 0)
			throw "Error size";
		pol_id = ++IDS;
		this->size = S;
		this->coef = new double[S];
		for (int i = 0; i < S; i++)
			coef[i] = C[i];
#ifdef _DEBUG
		std::cout << "The polynomial " << pol_id << " created" << std::endl;
#endif
	}

	polynomial::~polynomial()
	{
		size = 0;
		delete[] coef;
#ifdef _DEBUG
		std::cout << "The polynomial " << pol_id << " deleted" << std::endl;
#endif
	}

	polynomial::polynomial(const polynomial& other)
	{
		pol_id = ++IDS;
		this->size = other.size;
		if (other.size != 0)
			this->coef = new double[other.size];
		else
			this->coef = new double[1];
		for (int i = 0; i < other.size; i++)
			this->coef[i] = other.coef[i];
#ifdef _DEBUG
		std::cout << "Polynomial " << pol_id << " copied" << std::endl;
#endif
	}

	polynomial::polynomial(polynomial&& other) noexcept
	{
		std::swap(size, other.size);
		std::swap(coef, other.coef);
		pol_id = ++IDS;
#ifdef _DEBUG
		std::cout << "Вызван конструктор перемещения, объект с номером " << other.pol_id << " перемещен в объект " << pol_id << std::endl;
#endif
	}

	polynomial& polynomial::operator=(const polynomial& other)
	{
		pol_id = ++IDS;
		this->size = other.size;
		if (other.size != 0)
			this->coef = new double[other.size];
		else
			this->coef = new double[1];
		for (int i = 0; i < other.size; i++)
			this->coef[i] = other.coef[i];
#ifdef _DEBUG
		std::cout << "Polynomial " << pol_id << " copied" << std::endl;
#endif
		return *this;
	}

	polynomial& polynomial::operator=(polynomial&& other) noexcept
	{
		std::swap(size, other.size);
		std::swap(coef, other.coef);
		pol_id = ++IDS;
#ifdef _DEBUG
		std::cout << "Вызван конструктор перемещения, объект с номером " << other.pol_id << " перемещен в объект " << pol_id << std::endl;
#endif
		return *this;
	}

	polynomial& polynomial::operator+=(const polynomial& other)
	{
		if (this->size <= other.size)
		{
			double* x = new double[other.size];
			for (int i = 0; i < this->size; i++)
				x[i] = coef[i];
			for (int i = this->size; i < other.size; i++)
				x[i] = 0;
			delete[]coef;
			coef = new double[other.size];
			for (int i = 0; i < other.size; i++)
				coef[i] = x[i] + other.coef[i];
			delete[]x;
			this->size = other.size;
		}
		else
		{
			for (int i = 0; i < this->size; i++)
				coef[i] += other.coef[i];
		}
		return *this;
	}

	polynomial& polynomial::operator-=(const polynomial& other)
	{
		if (this->size <= other.size)
		{
			double* temp = new double[other.size];
			for (int i = 0; i < this->size; i++)
				temp[i] = coef[i];
			for (int i = this->size; i < other.size; i++)
				temp[i] = 0;
			delete[]coef;
			coef = new double[other.size];
			for (int i = 0; i < other.size; i++)
				coef[i] = temp[i] - other.coef[i];
			delete[]temp;
			this->size = other.size;
		}
		else
		{
			for (int i = 0; i < this->size; i++)
				coef[i] -= other.coef[i];
		}
		return *this;
	}

	polynomial& polynomial::operator*=(const polynomial& other)
	{
		double* temp = new double[this->size + other.size]();
		for (int i = 0; i < this->size; i++)
			for (int j = 0; j < other.size; j++)
				temp[i + j] += this->coef[i] * other.coef[j];
		delete[]this->coef;
		this->coef = new double[this->size + other.size]();
		this->coef = temp;
		temp=nullptr;
		this->size += other.size;

		return *this;
	}

	polynomial& polynomial::operator/=(const polynomial& other)
	{
		if (other.size > size)
			throw "Error! Size of the divisor > size of the divisible";
		double *temp = new double[size - other.size + 1]();
		polynomial temp1(size);
		for (int i = 0; i < size - other.size + 1; i++)
		{
			temp[this->size - other.size - i] = this->coef[this->size -1 - i] / other.coef[other.size -1];
			for (int j = 0; j < other.size; j++)
				temp1.coef[j + size - other.size-i] = other.coef[j] * temp[this->size - other.size - i];
			for (int j = 0; j < other.size; j++)
				this->coef[size - 1 - j-i] -= temp1.coef[size - 1 - j-i];
			for (int j = 0; j < size; j++)
				temp1.coef[j] = 0;
		}
		this->size = size - other.size + 1;
		this->coef = temp;
		temp = nullptr;
		
		return *this;
	}

	double polynomial::operator()(double k)
	{
		double result = 0;
		for (int i = 0; i < size; i++)
		{
			result += coef[i] * powl(k, i);
		}
		return result;
	}

	polynomial operator+(const polynomial& a, const polynomial& b)
	{
		return polynomial(a) += b;
	}

	polynomial operator-(const polynomial& a, const polynomial& b)
	{
		return polynomial(a) -= b;
	}

	polynomial operator*(const polynomial& a, const polynomial& b)
	{
		return polynomial(a) *= b;
	}

	polynomial operator/(const polynomial& a, const polynomial& b)
	{
		return polynomial(a) /= b;
	}

	std::ostream& operator<<(std::ostream& out, polynomial& other)
	{
		if (other.size != 0)
		{
			int k = other.size-1, count = 0;
			while ((other.coef[k] == 0) && (k>0))
				k--;
			if (k != 0)
			{
				for (int i = k; i != -1; i--)
				{
					if (other.coef[i] != 0)
					{
						out << other.coef[i] << "X^" << i;
						if (i > 0)
							if (other.coef[i - 1] > 0)
								out << "+";
					}						
				}
				out << std::endl;
				return out;
			}
			else
			{
				out << "Polinomial polynomial is zero" << std::endl;
				return out;
			}
		}
		else
		{
			out << "Polinomial is empty!" << std::endl;
			return out;
		}		
	}

}
