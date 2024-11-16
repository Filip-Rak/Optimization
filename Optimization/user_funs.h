#pragma once
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

matrix ff3T(matrix x, matrix ud1, matrix ud2);
