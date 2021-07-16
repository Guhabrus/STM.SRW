#pragma once
#include <fstream>
#include "Class.h"
class SensingElement;

void SolvingSystemOfAlgebraicEquations(SensingElement &ex1, SensingElement &ex2)
{
	
	std::cout << " Start Solving System Of Algebraic Equations " << std::endl;

	std::ofstream printFlux;		printFlux.open("ResultHeatFlux.txt");			printFlux << std::fixed;		printFlux.precision(3);			printFlux << "Q = [";
	std::ofstream printEmissivity;	printEmissivity.open("ResultEmissivity.txt");  printEmissivity << std::fixed;	printEmissivity.precision(3);	printEmissivity << " e = [";

	double a1, a2, b1, b2, c1, c2, x1, x2;
	double mistakeEmissivity = 0, mistakeHeatFlux = 0, MiddleHeatFlux = 0;
	
	for (int j = 1; j < SensingElement::sizeTim; j++)
	{


		a1 = -1;	a2 = -1;

		b1 = G * pow(ex1.T[j][0], 4);	b2 = G * pow(ex2.T[j][0], 4);

		c1 = -ex1.AbsorbedHeatFlux[j];	c2 = -ex2.AbsorbedHeatFlux[j];

		x1 = (c1*b2 - c2 * b1) / (a1*b2 - a2 * b1);
		x2 = (a1*c2 - a2 * c1) / (a1*b2 - a2 * b1);

		printFlux << x1 << ", ";
		printEmissivity << x2 << ", ";
		if (j >= 3)
		{
			mistakeEmissivity	+= fabs(Emissivity - x2);
			mistakeHeatFlux		+= fabs((ex1.GeneralHeatFlux[j] + ex2.GeneralHeatFlux[j]) / 2.0 - x1);
			MiddleHeatFlux		+= (ex1.GeneralHeatFlux[j] + ex2.GeneralHeatFlux[j]) / 2.0;
		}
		
	}

	
	printFlux		<< "];\n Time = [";
	printEmissivity << "];\n Time = [";

	for (int j = 1; j < SensingElement::sizeTim; j++)
	{
		printFlux << j * ex1.stepTime << ", ";
		printEmissivity << j * ex1.stepTime << ", ";
	}

	printFlux		<< "];\n plot(Time, Q,'b');\n hold on;\n grid on;\n";
	printEmissivity << "];\n plot(Time, e,'b');\n hold on;\n grid on;\n";



	std::cout << " Mistake of emissivity = " << mistakeEmissivity / (ex1.sizeTime - 3) << " \t Mistake of Heat Flux = " << mistakeHeatFlux / (ex1.sizeTime - 3) << std::endl;
	std::cout << " \nEnd Solving System Of Algebraic Equations " << std::endl;

	mistakeEmissivity = 100.0*(mistakeEmissivity / (ex1.sizeTime - 3))/Emissivity;		mistakeHeatFlux = 100.0*(mistakeHeatFlux / (ex1.sizeTime-3))/MiddleHeatFlux;

	std::cout << " Mistake of emissivity = " << mistakeEmissivity << " \t Mistake of Heat Flux = " << mistakeHeatFlux << std::endl;
	std::cout << " \nEnd Solving System Of Algebraic Equations " << std::endl;

}






#include <ctime>

double Shum(double T)
{

	//srand(static_cast<unsigned int>(time(NULL)));
	double random	= (rand() % 6 -3);
	
	double sigma	= del * T / 3;
	double f		= 1 / (sigma*sqrt(2 * M_PI))*exp(-(random*random / (2 * sigma*sigma)));

	if (f >= Min)
		return random;
	else
		return 0;
}


