#pragma once
#include <iostream>

namespace mathFraction
{
	class fraction
	{
	private:
		static unsigned IDS;
		unsigned fra_id;
		int integer;
		int denominator;
		int numerator;
	public:
		fraction(int integer, int numer, int denom);
		fraction(int numer, int denom);
		fraction(int integer);
		fraction();
		~fraction();
		fraction(const fraction& other);
		fraction(fraction&& other) noexcept;
		fraction& operator=(const fraction& other);
		fraction& operator=(fraction&& other) noexcept;
		fraction& operator+=(const fraction& other);
		fraction& operator-=(const fraction& other);
		fraction& operator*=(const fraction& other);
		fraction& operator/=(const fraction& other);
		friend fraction operator+(const fraction& a, const fraction& b);
		friend fraction operator-(const fraction& a, const fraction& b);
		friend fraction operator*(const fraction& a, const fraction& b);
		friend fraction operator/(const fraction& a, const fraction& b);
		friend std::ostream& operator<< (std::ostream& out, const fraction& fraction);
		int GreatestCommonDivisor(int numer, int denom);
		int LowestCommonDenominator(int a, int b);
		fraction& Reduction(int integer, int numer, int denom);
	};
}
