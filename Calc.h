#pragma once

#include "ConstData.h"
#include "Class.h"
#include "Functions.h"
#include "Input_Data.h"
#include "SolveInverse.h"




void Solve()
{

	SensingElement TZMK(0.05, 1200.0, 144.0);	TZMK.gran = 3000.0; TZMK.ac = 0.0000005;				/* Кф теплопроводности, Удельная теплоёмкость, плотность */
	TZMK.Print_InversHeatTransferProblemGradiend = true;

	TZMK.SetStartHeatFlux(1000.0);	TZMK.SetStartTemp(300.0);

	TZMK.SetNameOfSensingElement("TZMK");		TZMK.PrintCoefTemp();									/*Вычислениеи шагов вычислительно сетки*/
	
	double *X	= new double[n08];																		/*Подготовка данных для постоения сплайна*/
	double *Y	= new double[n08];
	double *Tx1 = new double[sT];

	for (int j = 0; j < n08; j++)
	{
		Y[j] = Tw08[j];
		X[j] = time08[j];
	}
	cubic_spline Example1;		Example1.build_spline(X, Y, n08);										 /*Построение сплайна внутри класса cubic_spline*/
	for (int j = 0; j < sT; j++)
		Tx1[j] = Example1.GetSplineNode(j*TZMK.GetStepTime());											/*Вывод точек сплайна в массив граничного условия*/




	int Board = 0;
	TZMK.DirectHeatTransferProblemFirstGY(Board, Tx1);	TZMK.PrintTempInX1Point();						/*Решение прямой задачи - определение поля температур*/
	
	

	TZMK.InversHeatTransferProblemGradiend();
	//TZMK.InversHeatTransferProblemFletcher();

	TZMK.PrintHeatFlux();
	TZMK.PrintParameters();
	TZMK.PrintGeneralHeatFlux();
	

	SensingElement KKK(2.2, 1105.0, 1900.0); KKK.gran = 3500.0;	KKK.ac = 0.00001;

	KKK.SetStartHeatFlux(30000.0);	KKK.SetStartTemp(320.0);

	KKK.Print_InversHeatTransferProblemGradiend = true;


	KKK.SetNameOfSensingElement("KKK");

	KKK = TZMK;



	KKK.DirectHeatTransferProblemFourthGY();	KKK.PrintTempInX1Point();

	KKK.InversHeatTransferProblemGradiend();

	KKK.PrintParameters();
	KKK.PrintHeatFlux();
	KKK.PrintGeneralHeatFlux();


	SolvingSystemOfAlgebraicEquations(TZMK, KKK);

	KKK.PrintInfo();
	
	/*

	

	
	
	
	KKK.DirectHeatTransferProblemSecondGY();
	KKK.PrintParameters();
	KKK.InversHeatTransferProblem();
	//KKK.InversHeatTransferProblemGradiend();

	SolvingSystemOfAlgebraicEquations(TZMK, KKK);*/

	delete[] X;
	delete[] Y;

}