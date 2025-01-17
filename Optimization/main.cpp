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
		// lab2();
		// lab3();
		// lab4();
		// lab5();
		lab6();
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

	// ---------- Table 1 and Table 2 ----------

	std::cout << "Solving for Table 1 and Table 2...";

	// Init random number generator
	srand(time(NULL));

	// Range for defined local minimums
	double x0_min = 1.0, x1_min = 1.0;  // Boundries for random number generation
	
	double sif_c = 0.05; // c size for internal 
	double sif_dc = 0.5; // c -> 0

	double sef_c = 5.0; // c size for external 
	double sef_dc = 2.0; // c -> inf

	// Common arguments
	double epsilon = 1e-3;
	int Nmax = 1000;

	//initial values for sym_NM
	matrix init_v_sym_NM = matrix(6, 1);

	init_v_sym_NM(0) = 1.0; //double side_size = 0.5;
	init_v_sym_NM(1) = 1.0; //double reflection_fator = 1.0;
	init_v_sym_NM(2) = 0.5; //double narrowing_factor = 0.5;
	init_v_sym_NM(3) = 2.0; //double expansion_factor = 2.0;
	init_v_sym_NM(4) = 0.5; //double reduction_factor = 0.5;
	init_v_sym_NM(5) = epsilon;
	const char delimiter = '\t';
	const string OUTPUT_PATH = "Output/lab_3/";
	ofstream tfun_file(OUTPUT_PATH + "out_1_tfun.txt");
	if (!tfun_file.is_open())
		return;

	const int k = 3;
	matrix value_a[k] = {matrix(1,1,4.0),matrix(1,1,4.4934),matrix(1,1,5)};
	
	for (int i = 0; i < k; i++) {
		for (int o = 0; o < 100; o++) {
			matrix x0(2, 1);
			x0(0) = x0_min + static_cast<double>(rand()) / RAND_MAX * (m2d(value_a[i]) - x0_min);
			x0(1) = x1_min + static_cast<double>(rand()) / RAND_MAX * (m2d(value_a[i]) - x1_min);
			double r_init = sqrt(pow(x0(0), 2) + pow(x0(1), 2));
			if (r_init > value_a[i])
			{
				o--;
				continue;
			}

			tfun_file << x0(0) << delimiter << x0(1) << delimiter;
			solution k = pen(SEF<0, ff3T>, x0, sef_c, sef_dc, epsilon, Nmax, init_v_sym_NM, value_a[i]);
			k.fit_fun(ff3T);
			tfun_file << k.x(0) << delimiter << k.x(1) << delimiter << sqrt(pow(k.x(0), 2) + pow(k.x(1), 2)) << delimiter;
			tfun_file << m2d(k.y) << delimiter << solution::f_calls << delimiter;
			solution::clear_calls();

			k = pen(SIF<1,ff3T>, x0, sif_c, sif_dc, epsilon, Nmax, init_v_sym_NM, value_a[i]);
			k.fit_fun(ff3T);
			tfun_file << k.x(0) << delimiter << k.x(1) << delimiter << sqrt(pow(k.x(0), 2) + pow(k.x(1), 2)) << delimiter;
			tfun_file << m2d(k.y) << delimiter << solution::f_calls << std::endl;

			if (k.flag == -2) std::cout << "-2!\n" << k.x(0) << " " << k.x(1) << "\n";
			solution::clear_calls();

		}
	}
	return;

	// ----- Real Problem ----- //
	std::cout << "\nSolving for real problem...";
	
	// Files
	std::ofstream opt_file(OUTPUT_PATH + "out_2_1_opt.txt");
	std::ofstream sim_file(OUTPUT_PATH + "out_2_2_sim.txt");

	// ----- Tab 3 -----

	// Starting conditions 
	double start_velocity = 5.0;	// [m/s]
	double omega = 10.0;	// [rad/s]
	matrix x0(2, new double[2]{start_velocity, omega});

	// Optimization
	double penalty_start = 100.0;
	double penalty_adjustment = 0.1;
	double penalty_multiplier = 1e7;
	solution opt_res = pen(ff3R, x0, penalty_start, penalty_adjustment, epsilon, Nmax, init_v_sym_NM, penalty_multiplier);

	// Save to file
	opt_file << opt_res << "\n";

	// ----- Simulation -----

	// Save optimized parameters
	double opt_velocity = opt_res.x(0);
	double opt_rotation = opt_res.x(1);

	// Starting parameters
	const double x_0 = 0;	// Starting horizontal position
	const double vx_0 = opt_res.x(0); // Starting horizontal speed
	const double y_0 = 100;	// Starting vertical position
	const double vy_0 = 0;	// Starting vertical speed
	const double omega_0 = opt_res.x(1); // Starting rotation
	matrix Y_0(4, new double[4] {x_0, vx_0, y_0, vy_0});

	// Time
	const double start_time = 0.0;
	const double end_time = 7.0;
	const double time_step = 0.01;

	// Resolve differential equation
	matrix* result = solve_ode(df3, start_time, time_step, end_time, Y_0, omega_0, NULL);

	// Save the result
	for (int i = 0; i < get_len(result[0]); ++i)
	{
		double t_i = result[0](i, 0);
		double X_i = result[1](i, 0);
		double Y_i = result[1](i, 2);

		sim_file << t_i << "\t" << X_i << "\t" << Y_i << "\n";
	}

	// Close the files
	opt_file.close();
	sim_file.close();

	// End debug with new line
	std::cout << "\n";
}

void lab4()
{
	// Common arguments
	double epsilon = 1e-3;
	int Nmax = 1000;

	// File output 
	const char delimiter = '\t';
	const string OUTPUT_PATH = "Output/lab_4/";
	const string INPUT_PATH = "Input/lab_4/";

	// ---------- Table 1 and Table 2 ----------

	std::cout << "Solving for Table 1 and Table 2...\n";

	// Open file to save test values
	ofstream tfun_file1(OUTPUT_PATH + "out_1_tfun.txt");
	if (!tfun_file1.good()) return;

	// Init random number generator
	srand(time(NULL));

	// Range for testing on both axis
	double max_values[2]{ 10.0, 10.0 };
	double min_values[2]{ -10.0, -10.0 };

	// Set step size for testing
	double steps[3]{ .05, .12, .0 };
	{
		// Test 100 random points with defined number of steps
		for (int step_i = 0; step_i < 3; step_i++)
		{
			for (int i = 0; i < 100; i++)
			{
				matrix test = matrix(2, new double[2] {
					min_values[0] + static_cast<double>(rand()) / RAND_MAX * (max_values[0] - min_values[0]),
						min_values[1] + static_cast<double>(rand()) / RAND_MAX * (max_values[1] - min_values[1])
					});
				tfun_file1 << test(0) << delimiter << test(1) << delimiter;

				// Steepest descent metohod -> step size determines the convergence of the function
				solution SSD = SD(ff4T, gradff4T, test, steps[step_i], epsilon, Nmax, NAN, NAN);
				SSD.fit_fun(ff4T, NAN, NAN);
				tfun_file1 << SSD.x(0) << delimiter << SSD.x(1) << delimiter << SSD.y(0) << delimiter;
				tfun_file1 << SSD.f_calls << delimiter << SSD.g_calls << delimiter;
				solution::clear_calls();

				// Conjugate gradient method -> guarantee of convergence, no need to calculate the Hessian matrix
				solution SCG = CG(ff4T, gradff4T, test, steps[step_i], epsilon, Nmax, NAN, NAN);
				SCG.fit_fun(ff4T, NAN, NAN);
				tfun_file1 << SCG.x(0) << delimiter << SCG.x(1) << delimiter << SCG.y(0) << delimiter;
				tfun_file1 << SCG.f_calls << delimiter << SCG.g_calls << delimiter;
				solution::clear_calls();

				// Newton's method -> guarantee of convergence, the Hessian matrix must be provided as additional Hf function.
				solution SNewton = Newton(ff4T, gradff4T, Hff4T, test, steps[step_i], epsilon, Nmax, NAN, NAN);
				SNewton.fit_fun(ff4T, NAN, NAN);
				tfun_file1 << SNewton.x(0) << delimiter << SNewton.x(1) << delimiter << SNewton.y(0) << delimiter;
				tfun_file1 << SNewton.f_calls << delimiter << SNewton.g_calls << delimiter << SNewton.H_calls << delimiter;
				solution::clear_calls();

				tfun_file1 << std::endl;
			}
		}
		tfun_file1.close();
	}
	try{
		matrix test_x = matrix(2, new double[2] {
			min_values[0] + static_cast<double>(rand()) / RAND_MAX * (max_values[0] - min_values[0]),
			min_values[1] + static_cast<double>(rand()) / RAND_MAX * (max_values[1] - min_values[1])
		});
		
		try {
			ofstream SD_tfun(OUTPUT_PATH + "out_SD_tfun.txt");
			if (!SD_tfun.good())  throw("FILE NO GOOD");
			for (int step_i = 0; step_i < 3; step_i++)
			{
				solution SSD = SD(ff4T, gradff4T, test_x, steps[step_i], epsilon, Nmax, 0, NAN);
				int* size_sd = get_size(SSD.ud);
				for (int i = 0; i < size_sd[1]; i++)
				{
					SD_tfun << SSD.ud(0, i) << delimiter << SSD.ud(1, i) << std::endl;
				}
				SD_tfun << SSD.x(0) << delimiter << SSD.x(1) << std::endl;
				SD_tfun << std::endl;
				delete[] size_sd;
				solution::clear_calls();
			}
			SD_tfun.close();
		} 
		catch (string ex_info)
		{
			throw ("SD_tfun:\n" + ex_info);
		}

		try {
			ofstream CG_tfun(OUTPUT_PATH + "out_CG_tfun.txt");
			if (!CG_tfun.good()) throw("FILE NO GOOD");
			for (int step_i = 0; step_i < 3; step_i++)
			{
				solution SCG = CG(ff4T, gradff4T, test_x, steps[step_i], epsilon, Nmax, 0, NAN);
				int* size_sd = get_size(SCG.ud);
				for (int i = 0; i < size_sd[1]; i++)
				{
					CG_tfun << SCG.ud(0, i) << delimiter << SCG.ud(1, i) << std::endl;
				}
				CG_tfun << SCG.x(0) << delimiter << SCG.x(1) << std::endl;
				CG_tfun << std::endl;
				delete[] size_sd;
				solution::clear_calls();
			}
			CG_tfun.close();
		}
		catch (string ex_info)
		{
			throw ("CG_tfun:\n" + ex_info);
		}

		try {
			ofstream Newton_tfun(OUTPUT_PATH + "out_Newton_tfun.txt");
			if (!Newton_tfun.good()) throw("FILE NO GOOD");
			for (int step_i = 0; step_i < 3; step_i++)
			{
				solution SNewton = Newton(ff4T, gradff4T, Hff4T, test_x, steps[step_i], epsilon, Nmax, 0, NAN);
				int* size_sd = get_size(SNewton.ud);
				for (int i = 0; i < size_sd[1]; i++)
				{
					Newton_tfun << SNewton.ud(0, i) << delimiter << SNewton.ud(1, i) << std::endl;
				}
				
				Newton_tfun << SNewton.x(0) << delimiter << SNewton.x(1) << std::endl;
				Newton_tfun << std::endl;
				delete[] size_sd;
				solution::clear_calls();
			}
			Newton_tfun.close();
		}
		catch (string ex_info)
		{
			throw ("Newton_tfun:\n" + ex_info);
		}

	}
	catch (string ex_info)
	{
		throw ("Graph values:\n" + ex_info);
	}

	// ---------- Real Problem ---------- //
	std::cout << "Solving for real problem...";

	// Input files
	ifstream x_matrix_file(INPUT_PATH + "XData.txt");
	ifstream y_matrix_file(INPUT_PATH + "YData.txt");

	// Matrix sizes
	const int cols = 100;
	const int x_rows = 3;
	const int y_rows = 1;

	// Load data from files into matrices (all integers)
	matrix x_matrix(x_rows, cols);	// Data per col [a] = { x_theta, x_grade1, x_grade2 }
	matrix y_matrix(y_rows, cols);	// Data per col [a] = { has_passed }

	x_matrix_file >> x_matrix;
	y_matrix_file >> y_matrix;

	// Close the input files
	x_matrix_file.close();
	y_matrix_file.close();

	matrix theta = {3, new double[3] {0.1f, 0.1f, 0.1f} };

	// Debug print with test values
	// std::cout << "J: " << get_cost(theta, y_matrix, x_matrix) << "\n";
	// std::cout << "Gradient:\n" << get_gradient(x_matrix, y_matrix, theta) << "\n";

	// Solve for various step lengths
	int real_n_max = 2e5;
	double real_epsilon = 1e-7;
	double step_length[] = { 0.01, 0.001, 0.0001 };
	int iterations = sizeof(step_length) / sizeof(double);

	// File reference for data output
	ofstream output_file;

	for (int i = 0; i < iterations; i++)
	{
		// Call the optimization function
		std::cout << "\n--- Step = " << step_length[i] << " # ";
		solution::clear_calls();
		solution opt_sol = CG(get_cost, get_gradient, theta, step_length[i], real_epsilon, real_n_max, y_matrix, x_matrix);
		std::cout << "#";

		// Open output file
		std::string out_number = to_string(i + 1);
		output_file.open(OUTPUT_PATH + "out_2_" + out_number + ".txt");

		// Output the solution to file
		output_file << opt_sol << "\n";
		output_file << "Accuracy: " << get_accuracy(opt_sol.x, x_matrix, y_matrix, cols) << "\n";

		for (int j = 0; j < cols; j++)
		{
			output_file << x_matrix(1, j) << delimiter << x_matrix(2, j) << delimiter
				<< y_matrix(0, j) << delimiter << get_h_l4(opt_sol.x, x_matrix, j) << "\n";
		}


		// Close output file
		output_file.close();
	}

	std::cout << "\n";
}

void lab5()
{
	// Common arguments
	double epsilon = 1e-3;
	int Nmax = 1000;

	// File output 
	const char delimiter = '\t';
	const string OUTPUT_PATH = "Output/lab_5/";
	const string INPUT_PATH = "Input/lab_5/";

	// Set weights
	init_weights();

	// Init random number generator
	srand(time(NULL));

	// ---------- Table 1 ----------
	std::cout << "Solving for Table 1...\n";

	double max_values[2]{ 10.0, 10.0 };
	double min_values[2]{ -10.0, -10.0 };

	{
		ofstream tfun_file1(OUTPUT_PATH + "out_1_tfun.txt");
		if (!tfun_file1.good()) return;

		for (int i = 0; i < 101; i++)
		{
			matrix test = matrix(2, new double[2] {
				min_values[0] + static_cast<double>(rand()) / RAND_MAX * (max_values[0] - min_values[0]),
					min_values[1] + static_cast<double>(rand()) / RAND_MAX * (max_values[1] - min_values[1])
				});
			set_weight(i);
			tfun_file1 << test(0) << delimiter << test(1) << delimiter;
			solution is1 = Powell(ff5T3_1, test, epsilon, Nmax, NAN, NAN);
			tfun_file1 << is1.x(0) << delimiter << is1.x(1) << delimiter << ff5T1_1(is1.x, NAN, NAN) << delimiter << ff5T2_1(is1.x, NAN, NAN) << delimiter << solution::f_calls << delimiter;
			solution::clear_calls();
			solution is10 = Powell(ff5T3_10, test, epsilon, Nmax, NAN, NAN);
			tfun_file1 << is10.x(0) << delimiter << is10.x(1) << delimiter << ff5T1_10(is10.x, NAN, NAN) << delimiter << ff5T2_10(is10.x, NAN, NAN) << delimiter << solution::f_calls << delimiter;
			solution::clear_calls();
			solution is100 = Powell(ff5T3_100, test, epsilon, Nmax, NAN, NAN);
			tfun_file1 << is100.x(0) << delimiter << is100.x(1) << delimiter << ff5T1_100(is100.x, NAN, NAN) << delimiter << ff5T2_100(is100.x, NAN, NAN) << delimiter << solution::f_calls << delimiter;
			solution::clear_calls();
			tfun_file1 << std::endl;
		}

		tfun_file1.close();
		//return;
	}

	/* Real Problem */
	std::cout << "Solving for Table 2...";

	// Output file
	ofstream rp_file(OUTPUT_PATH + "out2.txt");

	// Optimizer settings
	double real_epsilon = 1e-8;
	double real_nmax = 20000;

	// Constraints in meters
	// Beam length (l)
	const double beam_length_min = 0.200f;
	const double beam_length_max = 1.000f;

	// Cross sectional diameter (d)
	const double csd_min = 0.010;	
	const double csd_max = 0.050;

	// Loop over all weights
	const int weight_number = 101;
	for (int i = 0; i < weight_number; i++)
	{
		// Randomize beam's properties
		double beam_length = beam_length_min + static_cast<double>(rand()) / RAND_MAX * (beam_length_max - beam_length_min);
		double csd = csd_min + static_cast<double>(rand()) / RAND_MAX * (csd_max - csd_min);
		
		matrix x0 = matrix(2, new double[2] {beam_length, csd});

		// Set the weight
		set_weight(i);

		// Call the optimizer
		solution opt_res = Powell(ff5R, x0, real_epsilon, real_nmax, NULL, NULL);

		// Output data to the file
		rp_file
			<< get_weight() << delimiter			// Weight used for calculation
			<< beam_length * 1000.f << delimiter	// Random beam's length (l) in mm
			<< csd * 1000.f << delimiter;			// Random cross-sectional diameter (d) in mm

		rp_file 
			<< opt_res.x(0) * 1000.f << delimiter	// Optimized beam's length (l*) in mm
			<< opt_res.x(1) * 1000.f << delimiter	// Optimized cross-sectional diameter (d*) in mm
			<< ff5R_mass(opt_res.x)(0) << delimiter	// Mass for optimized X in kg
			<< ff5R_deflection(opt_res.x)(0) * 1000.f << delimiter	// Deflection for optimized X in mm
			<< opt_res.f_calls << delimiter			// Fit function calls to reach the solution
			<< opt_res.flag;						// Solution's flag

		// Add newline character between data
		rp_file << "\n";

		// Clear f_calls after saving to output
		solution::clear_calls();
	}

	// Close the output file
	rp_file.close();

	// Add new line for readability
	std::cout << "\n";
}
#include <memory>

class DoubleTabWrapper
{
public:
	double* ptr;
	DoubleTabWrapper(double* ptr)
	{
		this->ptr = ptr;
	}

	double* get() {
		return this->ptr;
	}
	
	~DoubleTabWrapper()
	{
		if(this->ptr)
			delete[] this->ptr;
	}
};

#include <typeinfo>
#include <thread>
void lab6()
{
	const string OUTPUT_PATH = "Output/lab_6/";
	const string INPUT_PATH = "Input/lab_6/";
	const char delimiter = '\t';
	matrix 
		lb(2, std::unique_ptr<double[]>(new double[2] { -5.0, -5.0 }).get()),
		ub(2, std::unique_ptr<double[]>(new double[2] { 5.0, 5.0 }).get());
	//return;
	double epsilon = 1e-2;
	int mi = 5;
	int lambd = 10;
	int Nmax = 1e+4;

	/* SKIP TO REAL PROBLEM */
	goto rp;

	//solution is1 = EA(ff6_T, 2, lb, ub, mi, lambd, 1, epsilon, Nmax);

	//goto tp;
	{
		ofstream tfun_file1(OUTPUT_PATH + "out_1_tfun.txt");
		if (!tfun_file1.good()) return;
		int pop = 0;
		double sigma_v[5] = { 0.01, 0.1, 1., 10., 100. };
		for (int i = 0; i < 5; i++)
		{
			//std::cout << sigma_v[i] << std::endl;
			for (int j = 0; j < 100; j++)
			{
				solution is1 = EA(ff6_T, 2, lb, ub, mi, lambd, sigma_v[i], epsilon, Nmax);
				tfun_file1 << is1.x(0) << delimiter
					<< is1.x(1) << delimiter 
					<< is1.y(0) << delimiter 
					<< solution::f_calls << delimiter 
					<< (is1.y(0) < epsilon) << std::endl;
				if (solution::f_calls > Nmax) pop++;
				solution::clear_calls();
				
			}
		}

		std::cout << "POP: " << pop << std::endl;
	}
tp:
	return;
	/* Real Problem */
rp:
	// Files
	ifstream ref_file(INPUT_PATH + "positions.txt");
	ofstream rp_opt_file(INPUT_PATH + "out_rp1_opt.txt");
	ofstream sim_file(OUTPUT_PATH + "out_rp2_sim.txt");

	// Optimizer settings
	double rp_lower_bound = 0.f;
	double rp_upper_bound = 3.f;
	matrix rp_lb(2, std::unique_ptr<double[]>(new double[2] {rp_lower_bound, rp_lower_bound}).get());
	matrix rp_ub(2, std::unique_ptr<double[]>(new double[2] {rp_upper_bound, rp_upper_bound}).get());

	double rp_epsilon = 1e-1, rp_sigma = 0.1f;
	int rp_mi = 5, rp_lambda = 10;
	int rp_n_max = 1e4, rp_n = 2;

	// Problem data
	int rows = 1001, cols = 2;
	matrix ref_data(rows, cols);
	ref_file >> ref_data;

	matrix rp_ud1 = NAN;
	matrix rp_ud2 = ref_data;

	// Optimize dampness values (b)
	solution::clear_calls();
	solution rp_solution = EA(ff6R, rp_n, rp_lb, rp_ub, rp_mi, rp_lambda, rp_sigma, rp_epsilon, rp_n_max, rp_ud1, rp_ud2);
	rp_opt_file << rp_solution;

	// Run simulation with optimized dampness
	matrix sim_res = ff6R_motion(rp_solution.x);
	sim_file << sim_res;

	// Print error of simulation
	matrix error = ff6R_error(sim_res, ref_data);
	std::cout << "RP_error: " << error << "\n";

	// Close files
	ref_file.close();
	rp_opt_file.close();
	sim_file.close();
}
