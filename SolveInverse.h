#include "Class.h"
#pragma once

class SensingElement;

typedef double funchion(double *, int);




inline void		SensingElement::DirectHeatTransferProblemFirstGY(int LeftBorder, double * T_GY1)
{
	std::cout << " Start Calc Direct Heat Transfer Problem First GY " << std::endl;
	double *P = new double[sizeX + 1];		//Прогоночный коэффициент альфа
	double *Q = new double[sizeX + 1];		//Прогоночный коэффициент бетта
	double a, b, c, d;
	double A = lamda / (Cp*density);

	for (int i = LeftBorder; i <= sizeX; i++)
		T[0][i] = StartTemp;

	for (int j = 1; j < sizeTime; j++)
	{
		for (int i = LeftBorder; i <= sizeX; i++)
		{
			if (i == LeftBorder)
			{
				P[i] = 0.0;
				Q[i] = T_GY1[j];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			}
			else
			{
				a = lamda / pow(stepX, 2);
				b = ((2 * lamda) / pow(stepX, 2)) + (Cp*density / stepTime);
				c = lamda / pow(stepX, 2);
				d = (-density * Cp*T[j - 1][i]) / stepTime;

				P[i] = a / (b - c * P[i - 1]);
				Q[i] = (c*Q[i - 1] - d) / (b - c * P[i - 1]);
			}
		}
		for (int i = sizeX; i >= LeftBorder; i--)
		{
			
			if (i == sizeX )
				T[j][i] = (2 * A*stepTime*lamda*Q[i - 1] - 2.*A*stepTime*stepX*HeatFluxRight + stepX * stepX*lamda*T[j - 1][i]) /
				(lamda*stepX*stepX + 2.*A*stepTime*lamda*(1 - P[i - 1]));
			else
				T[j][i] = P[i] * T[j][i + 1] + Q[i];
			/*
			if (i == sizeX)
				T[j][i] = 0.0;
			else
				T[j][i] = 0.0;*/
		}
	}


	std::ofstream Print;	Print.open("Tx1_True.txt");		Print << "Tx1 = [";

	for (int j = 0; j < sizeTime; j++) {
		Tx1True[j]	= T[j][X1] + Shum(T[j][X1]);
		Tx1[j]		= T[j][X1] + Shum(T[j][X1]);
	}
		


	for (int j = 0; j < sizeTime; j++)
	{

		Print << " Time = " << j * stepTime << " T = " << Tx1True[j] << " \n ";
	}
	Print.close();


	std::ofstream Print2;	Print2.open("T[0]_True.txt");	

	


	for (int j = 0; j < sizeTime; j++)
	{

		Print2 << " Time = " << j * stepTime << " T0 = " << Tx1True[j] << " \n ";
	}
	Print2.close();


	std::ofstream Print3;	Print3.open("Q_real.txt");	Print3 << "Q_real = [";

	for (int j = 1; j < sizeTime; j++)
		Print3 << (lamda * (T[j][0] - T[j][1]) / stepX + density * (stepX / 2.0)*(T[j][0] - T[j - 1][0]) / stepTime) + (pow(T[j][0],4)*G*Emissivity) << " , ";

	Print3 << " ];\n time = [ ";

	for (int j = 1; j < sizeTime; j++)
		Print3 << j * stepTime << " , ";

	Print3 << "];\n plot(time, Q_real,'r');\n hold on; \n grid on;";

	Print3.close();

	delete[] P;
	delete[] Q;

	std::cout << " End Calc Direct Heat Transfer Problem First GY " << std::endl;
}


inline void		SensingElement::DirectHeatTransferProblemSecondGY()
{
	double *P = new double[sizeX + 1];		//Прогоночный коэффициент альфа
	double *Q = new double[sizeX + 1];		//Прогоночный коэффициент бетта
	double a, b, c, d;	double Finish;
	double A = lamda / (Cp*density);

	for (int i = 0; i <= sizeX; i++)
		T[0][i] = StartTemp;
	int ptr = 0;
	std::cout << " 11111dc\n\n";
	Tx1True[0] = TempStart;
	for (int j = 1; j < sizeTime; j++)
	{
		

			for (int i = 0; i <= sizeX; i++)
			{

				if (i == 0)
				{
					P[i] = 1.0;
					Q[i] = (stepX / lamda)*(GeneralHeatFlux[j]-Emissivity * G*pow(1, 4));
				}
				else
				{
					a = lamda / pow(stepX, 2);
					b = ((2.0 * lamda) / pow(stepX, 2)) + (Cp*density / stepTime);
					c = lamda / pow(stepX, 2);
					d = (-density * Cp*T[j - 1][i]) / stepTime;

					P[i] = a / (b - c * P[i - 1]);
					Q[i] = (c*Q[i - 1] - d) / (b - c * P[i - 1]);
				}
			}
			for (int i = sizeX; i >= 0; i--)
			{
				if (i == sizeX)
					T[j][i] = (2 * A*stepTime*lamda*Q[i - 1] - 2.*A*stepTime*stepX*HeatFluxRight + stepX * stepX*lamda*T[j - 1][i]) /
					(lamda*stepX*stepX + 2.*A*stepTime*lamda*(1 - P[i - 1]));
				else
					T[j][i] = P[i] * T[j][i + 1] + Q[i];
			}



		Tx1[j] = T[j][X1];
		Tx1True[j] = T[j][X1];
	}
	
	std::ofstream print;	print.open("HeatFluxIn.txt");	print << " Q_in = [";

	for (int j = 0; j < sizeTime; j++)
		print << pow(2000. - j * stepTime, 2) / 500. + 5.*j*stepTime << ", ";
	print << " ];\n time = [";
	for (int j = 0; j < sizeTime; j++)
		print << j*stepTime<< ", ";
	print << " ];\n plot(time, Q_in, 'r');\n hold on;\n "; 
		

		print << " T_in = [";

		for (int j = 0; j < sizeTime; j++)
			print << Tx1True[j] << ", ";
		print << " ];\n time = [";
		for (int j = 0; j < sizeTime; j++)
			print << j * stepTime << ", ";
		print << " ];\n plot(time, T_in, 'r');\n hold on;\n ";
		print.close();

}	//


inline void		SensingElement::DirectHeatTransferProblemFourthGY()
{
	std::cout << " Start Calc Direct Heat Transfer Problem with Fourth GY for " << std::endl;
	double *P = new double[sizeX + 1];		//Прогоночный коэффициент альфа
	double *Q = new double[sizeX + 1];		//Прогоночный коэффициент бетта
	double a, b, c, d;	double Finish;
	double A = lamda / (Cp*density);

	for (int i = 0; i <= sizeX; i++)
		T[0][i] = StartTemp;
	int ptr = 0;
	Tx1True[0] = StartTemp;
	
	for (int j = 1; j < sizeTime; j++)
	{
		do
		{
			ptr++;

			if (ptr == FirstRun) { Finish = T[j - 1][0]; }
			else { Finish = T[j][0]; }

			for (int i = 0; i <= sizeX; i++)
			{

				if (i == 0)
				{
					P[i] = 1.0;
					Q[i] = (stepX / lamda)*(GeneralHeatFlux[j] - Emissivity * G*pow(Finish, 4));
				}
				else
				{
					a = lamda / pow(stepX, 2);
					b = ((2.0 * lamda) / pow(stepX, 2)) + (Cp*density / stepTime);
					c = lamda / pow(stepX, 2);
					d = (-density * Cp*T[j - 1][i]) / stepTime;

					P[i] = a / (b - c * P[i - 1]);
					Q[i] = (c*Q[i - 1] - d) / (b - c * P[i - 1]);
				}
			}
			for (int i = sizeX; i >= 0; i--)
			{
				if (i == sizeX)
					T[j][i] = (2 * A*stepTime*lamda*Q[i - 1] - 2.*A*stepTime*stepX*HeatFluxRight + stepX * stepX*lamda*T[j - 1][i]) /
					(lamda*stepX*stepX + 2.*A*stepTime*lamda*(1 - P[i - 1]));
				else
					T[j][i] = P[i] * T[j][i + 1] + Q[i];
			}


		} while (fabs(Finish - T[j][0]) > 1.1);
		ptr = 0;

		Tx1[j] = T[j][X1] + Shum(Tx1[j]);
		Tx1True[j] = T[j][X1] + Shum(Tx1[j]);
	}
	std::cout << " End Calc Direct Heat Transfer Problem with Fourth GY for " << std::endl;

}


inline void		SensingElement::InversHeatTransferProblem()
{
	/*
	double a, b, c, d;
	double *P = new double[sizeTime];
	double *Q = new double[sizeTime];

	for (int i = X1-1 ; i >= 0; i--)
	{
		for (int j = 0; j < sizeTime; j++)
		{
			if (j == 0)
			{
				a = 0.0;
				b = 1.0;
				c = 0.0;
				d = T[0][i];

				P[j] = 0.0;
				Q[j] = d / b;
			}
			else if (j > 0 && j < sizeTime - 1)
			{
				a = 1.0 / (2.*stepTime);
				b = -(lamda / (density*Cp)) / (pow(stepX, 2));
				c = -1.0 / (2.0*stepTime);
				d = (lamda / (Cp*density))*(T[j][i + 2] - 2.0*T[j][i + 1]) / pow(stepX, 2);

				P[j] = -a / (b + c * P[j - 1]);
				Q[j] = (d - c * Q[j - 1]) / (b + c * P[j - 1]);
			}
			else if (j == sizeTime - 1)
			{
				a = 0.0;
				b = 1 - (lamda / (Cp*density))*stepTime / (stepX*stepX);
				c = -1.0;
				d = stepTime * ((lamda / (Cp*density))*(T[j][i + 2] - 2.0*T[j][i + 1]) / pow(stepX, 2));

				P[j] = 0.0;
				Q[j] = (d - c * Q[j - 1]) / (b + c * P[j - 1]);
			}
			else { std::cout << " Mistake in OZT \n "; }
		}

		for (int j = sizeTime - 1; j >= 0; j--)
		{
			if (j == sizeTime - 1) {
				T[j][i] = Q[j];
			}
			else if (j >= 0 && j < sizeTime - 1) {
				T[j][i] = P[j] * T[j + 1][i] + Q[j];
			}
		}
	}



	delete[] P;
	delete[] Q;*/
	double A = lamda / (density*Cp);
	double p = stepX * stepX / (A*stepTime);

	for (size_t i = X1; i > 0; i--)
	{
		for (size_t j = 0; j < sizeTime - 1; j++)
			T[j + 1][i - 1] = (2 + p)*T[j + 1][i] - p * T[j][i] - T[j + 1][i + 1];
	}




	for (size_t j = 1; j < sizeTime; j++)
	{
		AbsorbedHeatFlux[j] = lamda * (T[j][0] - T[j][1]) / stepX + Cp * density*stepX / 2.*(T[j][0] - T[j - 1][0]) / stepTime;
	}

	for (size_t j = 1; j < sizeTime; j++)
		GeneralHeatFlux[j] = AbsorbedHeatFlux[j] + Emissivity * G*pow(T[j][0], 4);

}



inline void		SensingElement::InversHeatTransferProblemGradiend()
{
	std::cout << " Start Calc Invers Heat Transfer Problem Gradiend " << std::endl;

	for (size_t j = 1; j < sizeTime; j++)
		AbsorbedHeatFlux[j] = StartHeatFlux;

	double dQ = 1.0;
	double step_Q ;
	double dJ_dQ, dJ_dQ_plus_dt;
	double *dJ = new double[sizeTime];
	double f1, f2;
	std::cout << std::fixed;		std::cout.precision(7);
	double a = 0.0, b = 3500.0;

	double res1, res2, res3, res4;		double step_Q_left =a , step_Q_right = b, step_Q_middle;
	
	double y, z, fy, fz, fc, L;
	do
	{



		f1 = DirectHeatTransferProblemSecondGyForInvers();
		for (int j = 1; j < sizeTime; j++)
		{

			dJ_dQ = DirectHeatTransferProblemSecondGyForInvers();
			AbsorbedHeatFlux[j] += dQ;



			dJ_dQ_plus_dt = DirectHeatTransferProblemSecondGyForInvers();
			AbsorbedHeatFlux[j] -= dQ;



			dJ[j] = (dJ_dQ_plus_dt - dJ_dQ) / dQ;


		}
		step_Q_middle = (step_Q_left + step_Q_right) / 2.0;
		while (step_Q_right - step_Q_left > 0.4)
		{
			L = step_Q_right - step_Q_left;

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q_middle;
			fc = DirectHeatTransferProblemSecondGyForInvers();

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * step_Q_middle;



			y = step_Q_left + L / 4.0;	z = step_Q_right - L / 4.0;	

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * y;

			fy = DirectHeatTransferProblemSecondGyForInvers();

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * y;




			for(int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * z;

			fz = DirectHeatTransferProblemSecondGyForInvers();

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * z;


			if (fy < fc) {
				step_Q_right = step_Q_middle;
				step_Q_middle = y;
			}
			else if (fz < fc) {
				step_Q_left = step_Q_middle;
				step_Q_middle = z;
			}
			else if (fz > fc) {
				step_Q_left = y;
				step_Q_right = z;
			}


		}

		step_Q = (step_Q_left + step_Q_right) / 2.0;
		for (int j = 1; j < sizeTime; j++)
			AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q;

		f2 = DirectHeatTransferProblemSecondGyForInvers();


		step_Q_left = a;
		step_Q_right = b;

		std::cout << " J = " << fabs(f2 - f1) << " f2 = " << f2 << " Step_Q = " << step_Q << std::endl;
		
	} while (fabs(f2-f1) >ac);

	for (int j = 1; j < sizeTime; j++)
		GeneralHeatFlux[j] = AbsorbedHeatFlux[j] + Emissivity * G*pow(T[j][0], 4);






	if (Print_InversHeatTransferProblemGradiend)
	{
		std::ofstream Print;	Print.open("Tx_1_false.txt");	Print << "Tx_1_from_Inverse_Task=[";
		for (int j = 0; j < sizeTime; j++)
			Print << " Time = " << j * stepTime << "T = " << T[j][X1] << " \n ";

		Print << " Q = [";
		for (size_t j = 1; j < sizeTime; j++)
			Print << AbsorbedHeatFlux[j] << " , ";
		Print << "];\n time = [";


		for (size_t j = 1; j < sizeTime; j++)
			Print << j * stepTime << " , ";
		Print << "];\n plot(time, Q, 'b');\n hold on;\n grid on;";
	}
	
	delete[] dJ;

	std::cout << " End Calc Invers Heat Transfer Problem Gradiend " << std::endl;
}

inline void		SensingElement::InversHeatTransferProblemFletcher()
{
	std::cout << " Start Calc Invers Heat Transfer Problem Fletcher metod's " << std::endl;

	for (size_t j = 1; j < sizeTime; j++)
		AbsorbedHeatFlux[j] = StartHeatFlux;

	double dQ = 1.0;
	double step_Q;
	double dJ_dQ, dJ_dQ_plus_dt;

	double *dJ = new double[sizeTime];	double *d = new double[sizeTime];

	double f1, f2;
	std::cout << std::fixed;		std::cout.precision(5);
	double a = 0.0, b = 500.0;

	double res1, res2, res3, res4;		double step_Q_left = a, step_Q_right = b, step_Q_middle;

	double y, z, fy, fz, fc, L;

	int count = 1;	double dj_1 = 0, dj_2 = 0;	double betta=0;

	do {

		f1 = DirectHeatTransferProblemSecondGyForInvers();

		if (count == FirstRun) {
			
			for (int j = 1; j < sizeTime; j++)
			{

				dJ_dQ = DirectHeatTransferProblemSecondGyForInvers();
				AbsorbedHeatFlux[j] += dQ;



				dJ_dQ_plus_dt = DirectHeatTransferProblemSecondGyForInvers();
				AbsorbedHeatFlux[j] -= dQ;

				dJ[j] = (dJ_dQ_plus_dt - dJ_dQ) / dQ;

				d[j] = -dJ[j];

				
			}


			
			for (int j = 1; j < sizeTime; j++) {

				if (j == 1) {
					dj_1 = 1;
					dj_1 += dJ[j] * dJ[j];
				}
				else
					dj_1 += dJ[j] * dJ[j];
			
			}

			dj_1 = sqrt(dj_1);

			step_Q_middle = (step_Q_left + step_Q_right) / 2.0;

			while (step_Q_right - step_Q_left > 0.1) {

				L = step_Q_right - step_Q_left;

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * step_Q_middle;
				fc = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q_middle;



				y = step_Q_left + L / 4.0;	z = step_Q_right - L / 4.0;

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * y;

				fy = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - d[j] * y;




				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * z;

				fz = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - d[j] * z;


				if (fy < fc) {
					step_Q_right = step_Q_middle;
					step_Q_middle = y;
				}
				else if (fz < fc) {
					step_Q_left = step_Q_middle;
					step_Q_middle = z;
				}
				else if (fz > fc) {
					step_Q_left = y;
					step_Q_right = z;
				}


			}

			step_Q = (step_Q_left + step_Q_right) / 2.0;

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * step_Q;

			f2 = DirectHeatTransferProblemSecondGyForInvers();


			step_Q_left = a;
			step_Q_right = b;

			std::cout << " 1J = " << fabs(f2 - f1) << " f2 = " << f2 << " Step_Q = " << step_Q << std::endl;
			
		}
		else {
			for (int j = 1; j < sizeTime; j++)
			{

				dJ_dQ = DirectHeatTransferProblemSecondGyForInvers();
				AbsorbedHeatFlux[j] += dQ;


				dJ_dQ_plus_dt = DirectHeatTransferProblemSecondGyForInvers();
				AbsorbedHeatFlux[j] -= dQ;

				dJ[j] = (dJ_dQ_plus_dt - dJ_dQ) / dQ;

			}

			for (int j = 1; j < sizeTime; j++)
				dj_2 += dJ[j] * dJ[j];

			dj_2 = sqrt(dj_2);	betta = (dj_2*dj_2) / (dj_1*dj_1);	dj_1 = dj_2;

			for (int j = 1; j < sizeTime; j++)
				d[j] = -dJ[j] + betta * d[j];

			step_Q_middle = (step_Q_left + step_Q_right) / 2.0;

			while (step_Q_right - step_Q_left > 0.1) {
				L = step_Q_right - step_Q_left;

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * step_Q_middle;
				fc = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - d[j] * step_Q_middle;



				y = step_Q_left + L / 4.0;	z = step_Q_right - L / 4.0;

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * y;

				fy = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - d[j] * y;




				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * z;

				fz = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - d[j] * z;


				if (fy < fc) {
					step_Q_right = step_Q_middle;
					step_Q_middle = y;
				}
				else if (fz < fc) {
					step_Q_left = step_Q_middle;
					step_Q_middle = z;
				}
				else if (fz >= fc) {
					step_Q_left = y;
					step_Q_right = z;
				}


			}

			step_Q = (step_Q_left + step_Q_right) / 2.0;

			for (int j = 1; j < sizeTime; j++)
				AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + d[j] * step_Q;

			f2 = DirectHeatTransferProblemSecondGyForInvers();


			step_Q_left = a;
			step_Q_right = b;

			std::cout << " J = " << fabs(f2 - f1) << " f2 = " << f2 << " Step_Q = " << step_Q << std::endl;
			

			

		}
		count++;
	} while (fabs(f2 - f1) > AccuracuGradiend);

	for (int j = 1; j < sizeTime; j++)
		GeneralHeatFlux[j] = AbsorbedHeatFlux[j] + Emissivity * G*pow(T[j][0], 4);






	if (Print_InversHeatTransferProblemGradiend)
	{
		std::ofstream Print;	Print.open("Tx_1_false.txt");	Print << "Tx_1_from_Inverse_Task=[";
		for (int j = 0; j < sizeTime; j++)
			Print << " Time = " << j * stepTime << "T = " << T[j][X1] << " \n ";

		Print << " Q = [";
		for (size_t j = 1; j < sizeTime; j++)
			Print << AbsorbedHeatFlux[j] << " , ";
		Print << "];\n time = [";


		for (size_t j = 1; j < sizeTime; j++)
			Print << j * stepTime << " , ";
		Print << "];\n plot(time, Q, 'b');\n hold on;\n grid on;";
	}



	std::cout << " End Calc Invers Heat Transfer Problem Gradiend " << std::endl;

	delete[] dJ;	delete[] d;
	
}


inline double	SensingElement::DirectHeatTransferProblemSecondGyForInvers()
{
	double *P = new double[sizeX + 1];		//Прогоночный коэффициент альфа
	double *Q = new double[sizeX + 1];		//Прогоночный коэффициент бетта
	double a, b, c, d;	double Finish;
	double A = lamda / (Cp*density);

	for (size_t i = 0; i <= sizeX; i++)
		T[0][i] = StartTemp;

	for (int j = 1; j < sizeTime; j++)
	{
		for (int i = 0; i < sizeX; i++)
		{
			if (i == 0)
			{
				P[i] = 1;
				Q[i] = AbsorbedHeatFlux[j] * stepX / lamda;
			}
			else
			{
				a = lamda / pow(stepX, 2);
				b = ((2 * lamda) / pow(stepX, 2)) + (Cp*density / stepTime);
				c = lamda / pow(stepX, 2);
				d = (-density * Cp*T[j - 1][i]) / stepTime;

				P[i] = a / (b - c * P[i - 1]);
				Q[i] = (c*Q[i - 1] - d) / (b - c * P[i - 1]);
			}
		}
		for (int i = sizeX - 1; i >= 0; i--)
		{
			if (i == sizeX - 1)
				T[j][i] = (2 * A*stepTime*lamda*Q[i - 1] - 2.*A*stepTime*stepX*HeatFluxRight + stepX * stepX*lamda*T[j - 1][i]) /
				(lamda*stepX*stepX + 2.*A*stepTime*lamda*(1 - P[i - 1]));
			else
				T[j][i] = P[i] * T[j][i + 1] + Q[i];
		}
	}


	double J1 = 0;
	for (size_t j = 0; j < sizeTime; j++)
	{
		J1 += stepTime * pow((Tx1True[j] - T[j][X1]), 2);

		Tx1[j] = T[j][X1];
	}


	std::cout << std::fixed;
	std::cout.precision();
	//std::cout << "J = " << J1 << std::endl;

	delete[] P;
	delete[] Q;

	return J1;
}



inline SensingElement & SensingElement::operator=(const SensingElement & other)
{
	//if (this->GeneralHeatFlux != nullptr)
	//	delete[] this->GeneralHeatFlux;

	for (size_t j = 1; j < this->sizeTime; j++)
		this->GeneralHeatFlux[j] = other.GeneralHeatFlux[j];

	return *this;

}








inline void SensingElement::PrintParameters()
{
	std::ofstream print;	print.open(name + ".Temperatur field.txt");

	for (size_t j = 1; j < sizeTime; j++)
	{
		print << " Time = " << j * stepTime << std::endl;
		for (size_t i = 0; i <= sizeX; i++)
		{
			print << " | x = " << i * stepX << " T = " << T[j][i] << "| ";
		}
		print << "\n";
	}

	print.close();
}

inline void SensingElement::PrintHeatFlux() {
	std::ofstream print;	print.open(name + ".AbsorbedHeatFlux.txt");

	print << " Q = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << AbsorbedHeatFlux[j] << " , ";
	print << " ];\n time = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << j * stepTime << " , ";
	print << "];\n plot(time, Q,'b');\n hold on; \n grid on;\n";
	print.close();
}

inline void SensingElement::PrintTempInX1Point() {
	std::ofstream print;	print.open(name+".Температура в точке x1.txt");
	/*print << " Tx1true = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << Tx1True[j] << " , ";
	print << " ];\n time = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << j * stepTime << " , ";
	print << "];\n plot(time, Tx1true,'b');\n hold on; \n grid on;\n";
	*/
	print << " Tx1 = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << Tx1[j] << " , ";
	print << " ];\n time = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << j * stepTime << " , ";
	print << "];\n plot(time, Tx1,'b');\n hold on; \n grid on;\n";

	print.close();

}

inline void SensingElement::PrintGeneralHeatFlux()
{
	std::ofstream print;	print.open(name + ".GeneralHeatFlux.txt");

	print << " GenQ = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << GeneralHeatFlux[j] << " , ";
	print << " ];\n time = [";
	for (size_t j = 1; j < sizeTime; j++)
		print << j * stepTime << " , ";
	print << "];\n plot(time, GenQ,'b');\n hold on; \n grid on;\n";
	print.close();
}





/*

while (step_Q_right - step_Q_left > 0.01)
			{
				step_Q_middle = (step_Q_right + step_Q_left) / 2.0;

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q_right;
				res1 = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * step_Q_right;



				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q_left;
				res2 = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * step_Q_left;



				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * step_Q_middle;
				res3 = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * step_Q_middle;



				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] - dJ[j] * (step_Q_middle+1.0);
				res4 = DirectHeatTransferProblemSecondGyForInvers();

				for (int j = 1; j < sizeTime; j++)
					AbsorbedHeatFlux[j] = AbsorbedHeatFlux[j] + dJ[j] * (step_Q_middle + 1.0);

				if (((res1 - res2)*(res3 - res4)) < 0)
					step_Q_left = step_Q_middle;
				else
					step_Q_right = step_Q_middle;


			}*/