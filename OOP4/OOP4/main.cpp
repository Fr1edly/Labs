#include <iostream>
#include <math.h>
#include "fraction.h"
#include "polynomial.h"
#include <time.h>
#include "windows.h"
#include <iomanip>

using namespace mathFraction;
using namespace mathPolynom;

int main()
{
	srand(time(NULL));
	setlocale(LC_ALL, "ru");
	/*fraction a(-1,3);
	fraction b(-2,5);
	fraction c;
	std::cout << a;
	std::cout << b;
	c = a + b;
	std::cout << c;
	c = a - b;
	std::cout << c;
	c = a * b;
	std::cout << c;
	c = a / b;
	std::cout << c;*/
	double a1[]{ -36,-5,1 };
	polynomial a(3, a1);
	double b1[]{ -9,1 };
	polynomial b(2, b1);
	std::cout << a;
	std::cout << b;
	a /= b;
	std::cout << a;
}