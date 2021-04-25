#include <iostream>
#include <iomanip>
#include "fraction.h"
#include <string>
#include <cmath>


namespace mathFraction
{
	unsigned fraction::IDS = 0;

	fraction::fraction(int integer, int numer, int denom)
	{
		if (integer == 0 || numer == 0 || denom <= 0 || (integer < 0 && numer < 0))
			throw "Numerator or denominator size error!";

		fra_id = ++IDS;
		Reduction(integer, numer, denom);
#ifdef _DEBUG
		std::cout << "The fraction " << fra_id << " created" << std::endl;
#endif
	}

	fraction::fraction(int numer, int denom)
	{
		if (numer == 0 || denom <= 0)
			throw "Numerator or denominator size error!";

		int integer = 0;
		fra_id = ++IDS;
		Reduction(integer, numer, denom);
#ifdef _DEBUG
		std::cout << "The fraction " << fra_id << " created" << std::endl;
#endif
	}

	fraction::fraction(int integer)
	{
		if (integer == 0)
			throw "Size error!";
		this->integer = integer;
		this->numerator = 0;
		this->denominator = 0;
		fra_id = ++IDS;
#ifdef _DEBUG
		std::cout << "The fraction " << fra_id << " created" << std::endl;
#endif
	}

	fraction::fraction()
	{
		fra_id = ++IDS;
		int numer = rand() % 100 + 1;
		int denom = rand() % 100 + 1;
		int integer = 0;
		Reduction(integer, numer, denom);
#ifdef _DEBUG
		std::cout << "The fraction " << fra_id << " created" << std::endl;
#endif
	}

	fraction::~fraction()
	{
		this->integer = 0;
		this->numerator = 0;
		this->denominator = 0;
#ifdef _DEBUG
		std::cout << "The fraction " << fra_id << " deleted" << std::endl;
#endif
		fra_id = --IDS;
	}

	fraction::fraction(const fraction& other)
	{
		fra_id = ++IDS;
		this->integer = other.integer;
		this->numerator = other.numerator;
		this->denominator = other.denominator;
#ifdef _DEBUG
		std::cout << "Fraction " << fra_id << " copied" << std::endl;
#endif
	}

	fraction::fraction(fraction&& other) noexcept
	{
		fra_id = ++IDS;
		std::swap(integer, other.integer);
		std::swap(numerator, other.numerator);
		std::swap(denominator, other.denominator);
#ifdef _DEBUG
		std::cout << "Fraction " << other.fra_id << " moved to object " << fra_id << std::endl;
#endif
	}

	fraction& fraction::operator=(const fraction& other)
	{
		fra_id = ++IDS;
		this->integer = other.integer;
		this->numerator = other.numerator;
		this->denominator = other.denominator;
#ifdef _DEBUG
		std::cout << "Fraction " << fra_id << " copied" << std::endl;
#endif
		return *this;
	}

	fraction& fraction::operator=(fraction&& other) noexcept
	{
		fra_id = ++IDS;
		std::swap(integer, other.integer);
		std::swap(numerator, other.numerator);
		std::swap(denominator, other.denominator);
#ifdef _DEBUG
		std::cout << "Fraction " << other.fra_id << " moved to object " << fra_id << std::endl;
#endif
		return *this;
	}

	fraction& fraction::operator+=(const fraction& other)
	{
		int denom = other.denominator;
		int numer = other.numerator;
		if (denom == 0) denom = 1;
		if (this->denominator == 0) this->denominator = 1;
		int temp = LowestCommonDenominator(this->denominator, denom);
		int FactorForThis = temp / this->denominator;
		int FactorForOther = temp / denom;
		if (other.integer < 0)
			numer *= -1;
		if (this->integer < 0)
			this->numerator *= -1;
		this->numerator = (this->denominator * FactorForThis * this->integer + this->numerator * FactorForThis) + (denom * FactorForOther * other.integer + numer * FactorForOther);
		this->denominator *= FactorForThis;
		Reduction(0, this->numerator, this->denominator);
		return *this;
	}

	fraction& fraction::operator-=(const fraction& other)
	{
		int denom = other.denominator;
		int numer = other.numerator;
		if (denom == 0) denom = 1;
		if (this->denominator == 0) this->denominator = 1;
		int temp = LowestCommonDenominator(this->denominator, denom);
		int FactorForThis = temp / this->denominator;
		int FactorForOther = temp / denom;
		if (other.integer < 0)
			numer *= -1;
		if (this->integer < 0)
			this->numerator *= -1;
		this->numerator = (this->denominator * FactorForThis * this->integer + this->numerator * FactorForThis) - (denom * FactorForOther * other.integer + numer * FactorForOther);
		this->denominator *= FactorForThis;
		Reduction(0, this->numerator, this->denominator);
		return *this;
	}

	fraction& fraction::operator*=(const fraction& other)
	{
		int denom = other.denominator;
		int numer = other.numerator;
		if (denom == 0) denom = 1;
		if (this->denominator == 0) this->denominator = 1;
		if (other.integer < 0)
			numer *= -1;
		if (this->integer < 0)
			this->numerator *= -1;
		this->numerator = (this->integer * this->denominator + this->numerator) * (other.integer * denom + numer);
		this->denominator *= denom;
		Reduction(0, this->numerator, this->denominator);
		return *this;
	}

	fraction& fraction::operator/=(const fraction& other)
	{
		int denom = other.denominator;
		int numer = other.numerator;
		if (denom == 0) denom = 1;
		if (this->denominator == 0) this->denominator = 1;
		if (other.integer < 0)
			numer *= -1;
		if (this->integer < 0)
			this->numerator *= -1;
		this->numerator = (this->integer * this->denominator + this->numerator) * denom;
		this->denominator *= abs(other.integer * denom + numer);
		if (((numer < 0 || other.integer < 0) && (this->numerator<0 || this->integer<0)) || (numer < 0 || other.integer < 0) && (this->numerator > 0 || this->integer > 0))
			this->numerator *= -1;
		Reduction(0, this->numerator, this->denominator);
		return *this;
	}
	
	fraction operator+(const fraction& a, const fraction& b)
	{
		return fraction(a) += b;
	}

	fraction operator-(const fraction& a, const fraction& b)
	{
		return fraction(a) -= b;
	}

	fraction operator*(const fraction& a, const fraction& b)
	{
		return fraction(a) *= b;
	}

	fraction operator/(const fraction& a, const fraction& b)
	{
		return fraction(a) /= b;
	}

	std::ostream& operator<<(std::ostream& out, const fraction& fraction)
	{
		if (fraction.denominator == 0 && fraction.numerator == 0)
			out << fraction.integer << "\n";
		else
			if (fraction.integer == 0)
				out << fraction.numerator << "/" << fraction.denominator << "\n";
			else
				out << fraction.integer << ", " << fraction.numerator << "/" << fraction.denominator << "\n";
		return out;
	}

	int fraction::GreatestCommonDivisor(int numer, int denom)
	{
		numer = abs(numer);
		while (numer > 0 && denom > 0)
			if (numer > denom)
				numer %= denom;
			else
				denom %= numer;
		return numer + denom;
	}

	int fraction::LowestCommonDenominator(int a, int b)
	{
		return (a * b) / GreatestCommonDivisor(a, b);
	}

	fraction& fraction::Reduction(int integer, int numer, int denom)
	{
		int divisor = GreatestCommonDivisor(numer, denom);
		if (divisor > 1)
		{
			numer /= divisor;
			denom /= divisor;
		}
		if (abs(numer) >= denom)
		{
			if (numer % denom == 0)
			{
				integer += numer / denom;
				numer = 0;
				denom = 0;
			}
			else {
				integer += numer / denom;
				numer -= (numer / denom) * denom;
				if (integer < 0)
					numer = -numer;
			}
		}

		
		this->integer = integer;
		this->numerator = numer;
		this->denominator = denom;
		return *this;
	}
}
