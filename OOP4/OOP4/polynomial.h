#pragma once
#include <iostream>


namespace mathPolynom
{
	class polynomial
	{
	private:
		static unsigned IDS;
		unsigned pol_id;
		int size;
		double* coef;
	public:
		polynomial();
		polynomial(int S);
		polynomial(int S, double* C);
		~polynomial();
		polynomial(const polynomial&);
		polynomial(polynomial&&) noexcept;
		polynomial& operator=(const polynomial&);
		polynomial& operator=(polynomial&&) noexcept;
		polynomial& operator+=(const polynomial&);
		polynomial& operator-=(const polynomial&);
		polynomial& operator*=(const polynomial&);
		polynomial& operator/=(const polynomial&);
		friend polynomial operator+(const polynomial&, const polynomial&);
		friend polynomial operator-(const polynomial&, const polynomial&);
		friend polynomial operator*(const polynomial&, const polynomial&);
		friend polynomial operator/(const polynomial&, const polynomial&);
		double operator()(double);
		friend std::ostream& operator<< (std::ostream&, polynomial&);
	};
}