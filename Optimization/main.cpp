/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include <fstream>
#include <time.h>
#include <cstdlib>
#include "opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		// lab0();
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	// Shared arguments
	double epsilon = 1e-05;
	int n_max = 1000;

	// Expansion-specific arguments
	int x_min = -100, x_max = 100;  // Boundries for random number generation
	double d = 2, alpha = 2;

	// Lagrangea-specific arguments
	double gamma = 1e-200;

	// File output
	const char delimiter = ';';
	const string OUTPUT_PATH = "Output/";
	ofstream exp_file(OUTPUT_PATH + "out_exp.txt");
	ofstream fib_file(OUTPUT_PATH + "out_fib.txt");
	ofstream lag_file(OUTPUT_PATH + "out_lag.txt");

	// Init random number generator
	srand(time(NULL));

	// Table 1 and Table 2
	for (int i = 0; i < 100; i++)
	{
		// Narrow the range using expansion with random x0
		double x0 = x_min + (double)rand() / RAND_MAX * (x_max - x_min);
		double* range = expansion(ff1T, x0, d, alpha, n_max);

		if (exp_file.is_open())
			exp_file << x0 << delimiter << range[0] << delimiter << range[1] 
			<< delimiter << solution::f_calls << delimiter << "\n";
		
		// Use Fibonnaci's method
		solution::clear_calls();
		solution fib_result = fib(ff1T, range[0], range[1], epsilon);

		if (fib_file.is_open())
			fib_file << m2d(fib_result.x) << delimiter << m2d(fib_result.y) 
			<< delimiter << solution::f_calls << delimiter << fib_result.flag  << delimiter << "\n";

		// Use Lagrange's method
		solution::clear_calls();
		solution lag_result = lag(ff1T, range[0], range[1], epsilon, gamma, n_max);

		if (lag_file.is_open())
			lag_file << m2d(lag_result.x) << delimiter << m2d(lag_result.y) 
			<< delimiter << solution::f_calls << delimiter << lag_result.flag << delimiter<< "\n";


		// Deallocate memory
		delete[] range;
	}

	// Close the files
	exp_file.close();
	fib_file.close();
	lag_file.close();

	// Table 3
}

void lab1_1()
{
	// Shared arguments
	double epsilon = 0.001;
	int Nmax = 1000000;

	// Expansion Method
	double x0 = 30, d = 2, alpha = 2;
	double* p = expansion(ff1T, x0, d, alpha, Nmax);

	std::cout << " ( " << p[0] << " ; " << p[1] << " );\n";

	// Fibonacci's Method
	solution::clear_calls();

	std::cout << "-------- Metoda Fibonacciego --------" << "\n";
	std::cout << fib(ff1T, p[0], p[1], epsilon) << "\n";

	// Method based on Lagrange's Interpolation
	solution::clear_calls();
	double gamma = 1;

	std::cout << "-------- Metoda oparta na interpolacji Lagrange'a --------" << "\n";
	std::cout << lag(ff1T, p[0], p[1], epsilon, gamma, Nmax) << "\n";
	/*
		std::cout << "-------- P A I N --------" << "\n";

		solution fib_2 = fib(simulate_flow_temp, 1e-4, 1e-2, 1e-5);

		// va, vb, tb
		matrix Y0_1 = matrix(3, new double[3] {5, 1, 20});
		matrix* Yz1 = solve_ode(flow_and_temp, 0, 1, 2000, Y0_1, NULL, 0.005); // 0.005 = DA= 50CM
		std::cout << Yz1[1] << "\n";

		std::cout << "-------- P A I N 2 --------" << "\n";
		std::cout << "it fucking sucks: " << Yz1[1] << "\n";
	*/
	// solution opt = fib(simulate_flow_temp, -100, 100, epsilon, NULL, NULL);
	// std:cout << opt << "\n";

	delete[] p;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
