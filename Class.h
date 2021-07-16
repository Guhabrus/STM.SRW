#pragma once
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>


/*Начальная температура одинакова и равна 300 К. Задана руками в методах*/


#define FirstRun 1

//density - плотность 
class SensingElement
{
private:
	double lamda, Cp, density;
	
	double **T;

	double *Tx1, *Tx1True, *AbsorbedHeatFlux, *GeneralHeatFlux;

	std::size_t sizeX=sX, sizeTime = sT;

	double Length = 0.05, Time;

	double stepX ,	stepTime;

	double J = 0.;	double StartHeatFlux, StartTemp;
	
	std::string name;
	

public:
	bool Print_InversHeatTransferProblemGradiend;
	double gran, ac;

	void PrintParameters();
	void PrintHeatFlux();
	void PrintTempInX1Point();
	void PrintGeneralHeatFlux();
	void PrintInfo() {
		
		
			std::cout << std::fixed;	std::cout.precision(3);


			std::cout << " \t\t****Calc****\t\t\n";

			std::cout << " Step x = \t\t\t" << stepX << " \n step time = \t\t\t" << stepTime << "\n Time = \t\t\t" << Time << "\n x1 lenght \t\t= " << X1 * stepX << " x1 = \t\t\t " << X1 << std::endl;

			std::cout << "sizeX = \t\t\t" << sizeX << "\n ";
		

	}

	void PrintCoefTemp() { 
		double Fo = ((lamda / (density*Cp) * stepTime) / pow(((Length/sizeX)*X1), 2));

		std::cout << " Fo of ";	 printf("%s", &name[0]); std::cout << " = " << Fo << std::endl;
	}

	SensingElement(double lamda, double Cp, double density) {

		SetTime(TimeCalc);		SetStep();

		this->lamda = lamda;   this->Cp = Cp;   this->density = density;

		T = new double *[sizeTime];
		for (size_t i = 0; i < sizeTime; i++)
			T[i] = new double[sizeX+1];
	

		Tx1					= new double[sizeTime];
		Tx1True				= new double[sizeTime];
		AbsorbedHeatFlux	= new double[sizeTime];
		GeneralHeatFlux		= new double[sizeTime];

	};

	void SetStep() {
		stepX		= Length / (sizeX*1.0);
		stepTime	= Time / (sizeTime*1.0);
	}

	void SetNameOfSensingElement(std::string name) { this->name = name; }

	void SetStartHeatFlux(double StartHeatFlux) { this->StartHeatFlux = StartHeatFlux; }

	void SetStartTemp(double StartTemp) { this->StartTemp = StartTemp; }

	void DirectHeatTransferProblemFirstGY(int LeftBorder, double *T_GY1);
	

	void DirectHeatTransferProblemSecondGY();

	void DirectHeatTransferProblemFourthGY();

	void InversHeatTransferProblem();

	void InversHeatTransferProblemGradiend();
	void InversHeatTransferProblemFletcher();

	double  DirectHeatTransferProblemSecondGyForInvers();
	double	OptimazeStepGradiend(double step_Q);

	void SetSize(int sizeX, int sizeTime) { this->sizeX = sizeX;	this->sizeTime = sizeTime; }
	//void SetLengthOfTheSensingElement(double Length) { this->Length = Length; }
	void SetTime(double Time) { this->Time = Time; }
	
	double GetStepTime() { return this->stepTime; }

	friend void  SolvingSystemOfAlgebraicEquations(SensingElement &ex1, SensingElement &ex2);//Решение системы уравнений

	SensingElement &operator=(const SensingElement &other);


	~SensingElement() {

		PrintCoefTemp();

		for (int i = 0; i < sizeTime; i++)
			delete T[i];
		delete[] T;					T					= NULL;

		delete[] Tx1;				Tx1					= NULL;
		delete[] Tx1True;			Tx1True				= NULL;
		delete[] AbsorbedHeatFlux;	AbsorbedHeatFlux	= NULL;
		delete[] GeneralHeatFlux;	GeneralHeatFlux		= NULL;
	}

	static int sizeTim;

	//static int sizeX;
};
int SensingElement::sizeTim = sT;

//int SensingElement::sizeTime = sT;
//int SensingElement::sizeX = sX;




class Timer
{
public:
	Timer()
	{
		start = std::chrono::high_resolution_clock::now();
	}
	~Timer()
	{
		end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> difer = end - start;
		std::cout << " Duratioun - " << difer.count() << std::endl;
	}
private:
	std::chrono::time_point<std::chrono::steady_clock> start, end;
};

