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
		lab3();
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
	
	double sif_c = 5.0; // c size for internal 
	double sif_dc = 0.5; // c -> 0

	double sef_c = 5.0; // c size for external 
	double sef_dc = 2.0; // c -> inf

	// Common arguments
	double epsilon = 1e-3;
	int Nmax = 2000;

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
				if (x0(0) > x0(1))
					x0(0) -= r_init - m2d(value_a[i]);
				else
					x0(1) -= r_init - m2d(value_a[i]);
			}

			tfun_file << x0(0) << delimiter << x0(1) << delimiter;
			solution k = pen(SEF<0, ff3T >, x0, sef_c, sef_dc, epsilon, Nmax, init_v_sym_NM, value_a[i]);
			k.fit_fun(ff3T);
			tfun_file << k.x(0) << delimiter << k.x(1) << delimiter << sqrt(pow(k.x(0), 2) + pow(k.x(1), 2)) << delimiter;
			tfun_file << m2d(k.y) << delimiter << solution::f_calls << delimiter;
			solution::clear_calls();

			k = pen(SIF<0, ff3T>, x0, sif_c, sif_dc, epsilon, Nmax, init_v_sym_NM, value_a[i]);
			k.fit_fun(ff3T);
			tfun_file << k.x(0) << delimiter << k.x(1) << delimiter << sqrt(pow(k.x(0), 2) + pow(k.x(1), 2)) << delimiter;
			tfun_file << m2d(k.y) << delimiter << solution::f_calls << std::endl;
			solution::clear_calls();

		}
	}


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

}

void lab5()
{

}

void lab6()
{

}
