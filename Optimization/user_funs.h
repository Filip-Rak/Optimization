#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// Own Functions
// ------------------
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix flow_and_temp(double t, matrix Y, matrix ud1, matrix ud2);
matrix simulate_flow_temp(matrix x, matrix ud1, matrix ud2);