#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}


// Own Functions
// ------------------
matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	return pow(x, 2);
}

matrix flow_and_temp(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Result vector - changes in volume and temp
	matrix dY(3, 1);

	// Physical coefficients: viscosity, contraction, gravitational acceleration
	double a = 0.98, b = 0.63, g = 9.81;

	// A tank parameters
	double area_a = 0.5;
	double temp_a = 90.0;

	// B tank parameters
	double area_b = 1.0;
	double outflow_area_b = 0.00365665;

	// Inflow parameters
	double inflow_rate = 0.01;
	double temp_increase = 20.0;

	// Get outflows from each tank
	double outflow_a = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / area_a) : 0;
	double outflow_b = Y(1) > 1 ? a * b * outflow_area_b * sqrt(2 * g * Y(1) / area_b) : 0;

	// Changes in volume and temp:
	dY(0) = -outflow_a;	// Change in volume of tank A
	dY(1) = outflow_a + inflow_rate - outflow_b;	// Change in volume of tank B
	dY(2) = inflow_rate / Y(1) * (temp_increase - Y(2)) + outflow_a / Y(1) * (temp_a - Y(2));	// Change in temp in tank B

	return dY;
}

matrix simulate_flow_temp(matrix x, matrix ud1, matrix ud2)
{
	// Initial conditions
	double volume_a = 5, volume_b = 1, temp_b = 20;
	double conditions[3] = { volume_a, volume_b, temp_b };
	matrix Y0 = matrix(3, conditions);

	// Time
	double start_time = 0.0;
	double end_time = 2000.0;
	double time_step = 1.0;

	// Solve differential equation
	matrix* results = solve_ode(flow_and_temp, start_time, time_step, end_time, Y0, ud1, x);

	// Get max temp in B tank
	double max_temp_b = results[1](0, 2);

	for (int i = 1; i < get_len(results[0]); i++)
	{
		if (max_temp_b < results[1](i, 2))
			max_temp_b = results[1](i, 2);
	}

	// Get dabsolute deviation from the temp of 50
	double temp_deviation = abs(max_temp_b - 50);
	return temp_deviation;
}