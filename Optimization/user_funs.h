﻿#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include"ode_solver.h"


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// Own Functions LAB1
// ------------------

// f(x) = -cos(0.1x) * e^(-(0.1x-2*PI)^2) + 0.002 * (0.1x)^2
// @param matrix x - matrix 1x1 (scalar)
// @return matrix - matrix 1x1 (scalar)
matrix ff1T(matrix, matrix = NAN, matrix = NAN);

matrix flow_and_temp(double t, matrix Y, matrix ud1, matrix ud2);
matrix simulate_flow_temp(matrix x, matrix ud1, matrix ud2);

// Own Functions LAB2
// ------------------

// f(x) = x1^2 + x2^2 - cos(2.5*PI*x1) - cos(2.5*PI*x2) + 2
// @param matrix x - matrix 1x2 (vertical vector)
// @return matrix - matrix 1x1 (scalar)
matrix ff2T(matrix, matrix = NAN, matrix = NAN);

matrix df2(double t, matrix Y, matrix ud1, matrix ud2);

matrix ff2R(matrix x, matrix k1, matrix k2);

// Own Functions LAB3
// ------------------

// @param matrix x - matrix 1x2 (vertical vector)
matrix g1(matrix x, matrix = NAN);
// @param matrix x - matrix 1x2 (vertical vector)
matrix g2(matrix x, matrix = NAN);
// @param matrix x - matrix 1x2 (vertical vector)
matrix g3(matrix x, matrix a);


typedef matrix(*g_fun)(matrix, matrix);
constexpr int n_sets = 2;
constexpr int n_fun = 3;
/// @summary usage:
/// @pre set number of sets (n_sets)
///	@pre set number of functions in all steps (empty values must be set to nullptr!) 
///	@pre specify g functions and order 
///	@example 
/// n_sets = 2;
/// n_fun = 4;
/// gl_g_tab[n_sets][n_fun] = { { g1, nullptr, g3, g2 }, { nullptr, g2, nullptr, nullptr } }
constexpr g_fun gl_g_tab[n_sets][n_fun] = { { g1,g2,g3,}, {g3,g2,g1} };

/// usage: SZF<number of the set of g functions inside gl_g_tab, function for calculating F>
/// @param matrix x - point to calculate
/// @param matrix c_and_other - vertical vector (n,1), n >= 1, c_and_other(0) must be equal to c
template <int k, matrix(*f)(matrix, matrix, matrix) >
matrix SEF(matrix x, matrix c_and_other, matrix ud1 = NAN)
{
	if (k >= n_sets || k < 0)
		throw("matrix SEF(matrix x, matrix ud1, matrix ud2): nie ma takiego zbioru" + k);
	double sum = 0;
	for (int i = 0; i < n_fun; i++)
	{
		if (gl_g_tab[k][i] != nullptr)
		sum += pow(max(0.0,m2d(gl_g_tab[k][i](x, c_and_other(1)))),2);
	}
	//std::cout << "SEF: " << (f(x, NAN, NAN) + sum * c_and_other(0)) << endl;
	return (f(x, NAN, NAN) + sum * c_and_other(0));
}

static double max_value = 0.0;

/// usage: SWF<number of the set of g functions inside gl_g_tab, function for calculating F>
/// @param matrix x - point to calculate
/// @param matrix c_and_other - vertical vector (n,1), n >= 1, c_and_other(0) must be equal to c
template <int k, matrix(*f)(matrix, matrix, matrix)>
matrix SIF(matrix x, matrix c_and_other, matrix ud1 = NAN)
{
	if (k >= n_sets || k < 0)
		throw("matrix SIF(matrix x, matrix ud1, matrix ud2): nie ma takiego zbioru " + k);
	double sum = 0;
	for (int i = 0; i < n_fun; i++)
	{
		if(gl_g_tab[k][i] != nullptr)
		sum += 1.0/fabs(m2d(gl_g_tab[k][i](x, c_and_other(1))));
	}
	//std::cout << "SIF: " << sum * c_and_other(0) << endl;
	
	return (f(x, NAN, NAN) - sum * c_and_other(0));
}

matrix ff3T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix df3(double t, matrix Y, matrix ud1, matrix ud2);

matrix ff3R(matrix x, matrix ud1, matrix ud2);

// Own Functions LAB4
// ------------------


matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix gradff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix Hff4T(matrix x, matrix ud1, matrix ud2);

double get_h_l4(matrix theta, matrix X, int col);

matrix get_cost(matrix X, matrix Y, matrix theta);

matrix get_gradient(matrix X, matrix Y, matrix theta);

matrix get_accuracy(matrix theta, matrix X, matrix Y, int cols);

// Own Functions LAB5
// ------------------

void set_weight(int i);
void init_weights();
double get_weight();

/* Test Function */

matrix ff5T1_1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T2_1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T3_1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T1_10(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T2_10(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T3_10(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T1_100(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T2_100(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T3_100(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

/* Real Problem */

/*
	@brief Calculates beam's stress
	@param matrix x(0) = beam length (meters)
	@param matrix x(1) = cross-sectional diameter (meters)
	@return stress in Pa = N/m^2 as 1d matrix
*/
matrix ff5R_stress(matrix x);

/*
	@brief Calculates beam's mass
	@param matrix x(0) = beam length (meters)
	@param matrix x(1) = cross-sectional diameter (meters)
	@return mass in kg as 1d matrix
*/
matrix ff5R_mass(matrix x);

/*
	@brief Calculates beam's deflection
	@param matrix x(0) = beam length (meters)
	@param matrix x(1) = cross-sectional diameter (meters)
	@return deflection in meters as 1d matrix
*/
matrix ff5R_deflection(matrix x);

/*
	@brief Calculates penalty based on stress and deflection
	@param matrix x(0) = beam length (meters)
	@param matrix x(1) = cross-sectional diameter (meters)
	@return penalty in meters as 1d matrix
*/
matrix ff5R_penalty(matrix x);

/*
	@brief Beam model with applied penalty
	@param matrix x(0) = beam length (meters)
	@param matrix x(1) = cross-sectional diameter (meters)
	@return fitness/error scalar as 1d matrix
*/
matrix ff5R(matrix x, matrix ud1 = NULL, matrix ud2 = NULL);

// Own Functions LAB6
// ------------------

matrix ff6_T(matrix x, matrix ud1 = NULL, matrix ud2 = NULL);

/* Real Problem */
matrix ff6R_derivative(double t, matrix Y, matrix dampness, matrix ud2 = NAN);
matrix ff6R_motion(matrix x);
double ff6R_error(matrix simulation_results, matrix experimental_data);
matrix ff6R(matrix x, matrix ud1, matrix ud2);