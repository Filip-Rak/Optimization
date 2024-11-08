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
		// lab1();
		lab2();
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
	// Common arguments
	double epsilon = 1e-05;
	int n_max = 1000;

	// Expansion-specific arguments
	int x_min = -100, x_max = 100;  // Boundries for random number generation
	double d = 2, alpha = 13;

	// Lagrangea-specific arguments
	double gamma = 1e-200;

	// File output
	const char delimiter = '\t';
	const string OUTPUT_PATH = "Output/lab_1/";
	ofstream exp_file(OUTPUT_PATH + "out_1_1_exp.txt");
	ofstream fib_file(OUTPUT_PATH + "out_1_2_fib.txt");
	ofstream lag_file(OUTPUT_PATH + "out_1_3_lag.txt");


	// ---------- Table 1 and Table 2 ----------

	// Init random number generator
	srand(time(NULL));

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
			<< delimiter << solution::f_calls << delimiter << fib_result.flag << delimiter << "\n";

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

	// ---------- Graph ----------
	
	// Open new files for output
	fib_file.open(OUTPUT_PATH + "out_2_1_fib.txt");
	lag_file.open(OUTPUT_PATH + "out_2_2_lag.txt");

	// Fixed range
	double a = -100, b = 100;

	// Fibonacci function call
	solution::clear_calls();
	solution fib_result = fib(ff1T, a, b, epsilon);
	fib_file << m2d(fib_result.x) << delimiter << m2d(fib_result.y)
		<< delimiter << solution::f_calls << delimiter << fib_result.flag
		<< delimiter << "\n\n" << fib_result.ud << "\n";

	// Lagrange Functions call
	solution::clear_calls();
	solution lag_result = lag(ff1T, a, b, epsilon, gamma, n_max);
	lag_file << m2d(lag_result.x) << delimiter << m2d(lag_result.y)
		<< delimiter << solution::f_calls << delimiter << lag_result.flag
		<< delimiter << "\n\n" << lag_result.ud << "\n";

	// Close the files
	fib_file.close();
	lag_file.close();

	// ---------- Table 3 ----------

	// Open new files
	fib_file.open(OUTPUT_PATH + "out_3_1_fib.txt");
	lag_file.open(OUTPUT_PATH + "out_3_2_lag.txt");
	
	// Table specific function arguments
	double range[] = { 1e-4, 1e-2 };

	// Call Fibonacci
	solution::clear_calls();
	fib_result = fib(simulate_flow_temp, range[0], range[1], epsilon);
	fib_file << m2d(fib_result.x) << delimiter << m2d(fib_result.y) << delimiter
		<< solution::f_calls << delimiter << fib_result.flag << delimiter << "\n";

	// Call Lagrange
	solution::clear_calls();
	lag_result = lag(simulate_flow_temp, range[0], range[1], epsilon, gamma, n_max);
	lag_file << m2d(lag_result.x) << delimiter << m2d(lag_result.y) << delimiter 
		<< solution::f_calls << delimiter << lag_result.flag << delimiter << "\n";

	// Close the files
	fib_file.close();
	lag_file.close();

	// ---------- Simulation ----------
	// Initial conditions
	double volume_a = 5, volume_b = 1, temp_b = 20;
	double conditions[3] = { volume_a, volume_b, temp_b };
	matrix Y0 = matrix(3, conditions);

	// Time
	double start_time = 0.0;
	double end_time = 2000.0;
	double time_step = 1.0;

	// Solve differential equation with result from Fibonacci
	matrix* ode_fib_result = solve_ode(flow_and_temp, start_time, time_step, end_time, Y0, NULL, fib_result.x(0));
	
	// Solve differential equation with result from Lagrange
	matrix* ode_lag_result = solve_ode(flow_and_temp, start_time, time_step, end_time, Y0, NULL, lag_result.x(0));

	// Open files and save results
	fib_file.open(OUTPUT_PATH + "out_4_1_fib.txt");
	lag_file.open(OUTPUT_PATH + "out_4_2_lag.txt");

	fib_file << ode_fib_result[1];
	lag_file << ode_lag_result[1];

	// Close the files
	fib_file.close();
	lag_file.close();
}

void lab2()
{
	// Common arguments
	double epsilon = 1e-06;
	int n_max = 1000;
	double contr = 0.5, expa = 2; // for HJ: contr -> alpha, for Rosen: expa -> alpha, contr -> beta

	// ---------- Table 1 and Table 2 ----------
	
	std::cout << "Solving for Table 1 and Table 2...";

	// Init random number generator
	srand(time(NULL));

	// Range for defined local minimums
	double x_min = -1.0, x_max = 1.0;  // Boundries for random number generation

	// Set step size for testing
	double step[3] = { .75, .5, .25 }; // Step size for each iteration

	// File output
	const char delimiter = '\t';
	const string OUTPUT_PATH = "Output/lab_2/";
	ofstream tfun_file(OUTPUT_PATH + "out_1_tfun.txt");

	// Check if the file has been correctly opened
	if (!tfun_file.is_open())
		return;

	bool found_candidate = false;
	solution H0, R0;

	for (int j = 0; j < 3; j++)
	{

		for (int i = 0; i < 100; i++)
		{
			// Draw a 2D point
			matrix x(2,1);
			x(0) = x_min + static_cast<double>(rand()) / RAND_MAX * (x_max - x_min);
			x(1) = x_min + static_cast<double>(rand()) / RAND_MAX * (x_max - x_min);
			
			// Save points in the file
			tfun_file << x(0) << delimiter << x(1) << delimiter;
			
			// Calculate and write found answer to the same file, Hooke-Jeeves method.
			solution y0 = HJ(ff2T, x, step[j], contr, epsilon, n_max);
			tfun_file << y0.x(0) << delimiter << y0.x(1) << delimiter
				<< m2d(y0.y) << delimiter << solution::f_calls << delimiter 
				<< (abs(m2d(y0.y)) < epsilon) << delimiter;
			solution::clear_calls();

			// Calculate and write found answer to the same file, Rosenbrock method.
			solution y1 = Rosen(ff2T, x, matrix(2, 1, step[j]), expa, contr, epsilon, n_max);
			tfun_file << y1.x(0) << delimiter << y1.x(1) << delimiter
				<< m2d(y1.y) << delimiter << solution::f_calls << delimiter 
				<< (abs(m2d(y1.y)) < epsilon);
			solution::clear_calls();
			if (!found_candidate && (abs(m2d(y1.y)) < epsilon) && (abs(m2d(y0.y)) < epsilon))
			{
				found_candidate = true;
				H0 = y0;
				R0 = y1;
			}

			tfun_file << std::endl;
		}

	}
	// Close file
	tfun_file.close();

	ofstream graph_file(OUTPUT_PATH + "out_2_tfun.txt");
	if (!graph_file.is_open())
		return;
	int* hl = get_size(H0.ud);
	for (int i = 0; i < hl[1]; i++)
		graph_file << H0.ud(0, i) << delimiter << H0.ud(1, i) << std::endl;
	
	graph_file << "\n-----\n";
	delete[] hl;
	hl = get_size(R0.ud);
	for (int i = 0; i < hl[1]; i++)
		graph_file << R0.ud(0, i) << delimiter << R0.ud(1, i) << std::endl;
	
	delete[] hl;
	graph_file.close();
	return;
	// ---------- Table 3 ----------

	std::cout << "\nSolving for Table 3...";

	// Open files
	ofstream hook_file(OUTPUT_PATH + "out_3_1_hook.txt");
	ofstream rosen_file(OUTPUT_PATH + "out_3_2_rosen.txt");
	
	// Problem parameters
	double starting_step = 1.0;
	double k_values[2] = { 2.0, 6.0 };
	matrix x0 = matrix(2, k_values);

	// Save the results of Hook's method
	solution hook_opt = HJ(ff2R, x0, starting_step, contr, epsilon, n_max);
	hook_file 
		<< hook_opt.x(0) << delimiter  << hook_opt.x(1) << delimiter
		<< m2d(hook_opt.y) << delimiter << hook_opt.f_calls << delimiter 
		<< hook_opt.flag;
	solution::clear_calls();

	// Save the results of Rosenbrock's method
	solution rosen_opt = Rosen(ff2R, x0, matrix(2, 1, starting_step), expa, contr, epsilon, n_max);
	rosen_file
		<< rosen_opt.x(0) << delimiter << rosen_opt.x(1) << delimiter
		<< m2d(rosen_opt.y) << delimiter << rosen_opt.f_calls << delimiter
		<< rosen_opt.flag;
	solution::clear_calls();

	// Close the files
	hook_file.close();
	rosen_file.close();

	// ---------- Simulation ----------

	std::cout << "\nSolving for Simulation table...\n";

	// Open files
	hook_file.open(OUTPUT_PATH + "out_4_1_hook.txt");
	rosen_file.open(OUTPUT_PATH + "out_4_2_rosen.txt");

	// Initial condition
	matrix Y0 = matrix(2, new double[2] {0.0, 0.0});	// Start angle and speed

	// Time parameters
	double start_time = 0.0;
	double end_time = 100.0;
	double time_step = 0.1;

	// Solve for Hook's optimization
	matrix ud1 = hook_opt.x(0);	// K1
	matrix ud2 = hook_opt.x(1);	// K2
	matrix* hook_results = solve_ode(df2, start_time, time_step, end_time, Y0, ud1, ud2);
	hook_file << hook_results[1];

	// Solve for Rosen's optimization
	ud1 = rosen_opt.x(0); // K1
	ud2 = rosen_opt.x(1); // K2
	matrix* rosen_results = solve_ode(df2, start_time, time_step, end_time, Y0, ud1, ud2);
	rosen_file << rosen_results[1];

	// Close the files
	hook_file.close();
	rosen_file.close();
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
