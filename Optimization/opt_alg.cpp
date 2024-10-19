#include"opt_alg.h"
#include <vector>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };

		//Tu wpisz kod funkcji
		// --------------------
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		if (X0.y == X1.y)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}

		if (X1.y > X0.y)
		{
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d;
				return p;
			}
		}

		int i = 1;
		solution X2(X1);	// "x + i"
		while (true)
		{
			if (solution::f_calls > Nmax)
			{
				// "Error"
				throw("Nie znaleziono przedzialu\n");
			}
			
			X1 = X2;
			X2.x = X0.x + pow(alpha, i++) * d;
			X2.fit_fun(ff, ud1, ud2);
			// X1.fit_fun(ff, ud1, ud2);

			if (!(X1.y <= X2.y))
				break;
		}

		if (d > 0)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X2.x);
			return p;
		}

		p[0] = m2d(X2.x);
		p[1] = m2d(X0.x);
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		std::vector<int> fibs;
		fibs.push_back(1);
		fibs.push_back(1);

		int new_fib = fibs.back();
		while (new_fib <= (b - a) / epsilon)
		{
			new_fib = fibs[fibs.size() - 1] + fibs[fibs.size() - 2];
			fibs.push_back(new_fib);
		}

		int k = fibs.size() - 1;

		double a0 = a, b0 = b;
		solution c0(b0 - (double)fibs[k - 1] / fibs[k] * (b0 - a0));
		solution d0((double)fibs[k - 1] / fibs[k] * (b0 - a0));

		for (int i = 0; i <= k - 3; i++)
		{
			if (c0.fit_fun(ff, ud1, ud2) < d0.fit_fun(ff, ud1, ud2))
				b0 = m2d(d0.x);
			else
				a0 = m2d(c0.x);

			c0.x = b0 - (double)fibs[k - 1] / fibs[k] * (b0 - a0);
			d0.x = a0 + (double)fibs[k - i - 2] / fibs[k - i - 1] * (b0 - a0);
		}

		Xopt = c0;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (const std::string& ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double a0 = a, b0 = b, c0 = (a + b) / 2;
		double dm1 = 0;
		double d0 = 0;

		while (true)
		{
			double fa = m2d(ff(a0, ud1, ud2));
			double fb = m2d(ff(b0, ud1, ud2));
			double fc = m2d(ff(c0, ud1, ud2));

			double l = fa * (b0 * b0 - c0 * c0) + fb * (c0 * c0 - a0 * a0) + fc * (a0 * a0 - b0 * b0);
			double m = fa * (b0 - c0) + fb * (c0 - a0) + fc * (a0 - b0);

			if (m <= 0)
				throw ("Brak rozwiązania\n");

			dm1 = d0;
			d0 = 0.5 * l / m;
			if (a0 < d0 && d0 < c0)
			{
				if (ff(d0, ud1, ud2) < ff(c0, ud1, ud2))
				{
					c0 = d0;
					b0 = c0;
				}
				else
					a0 = d0;

			}
			else if (c0 < d0 && d0 < b0)
			{
				if (ff(d0, ud1, ud2) < ff(c0, ud1, ud2))
				{
					a0 = c0;
					c0 = d0;
				}
				else
					b0 = d0;
			}
			else
			{
				//throw ("Brak rozwiązania\n");
			}


			i++;
			if (i > Nmax)
				throw ("Nmax");

			if (b0 - a0 < epsilon || fabs(d0 - dm1) < gamma)
				break;

		}

		Xopt.x = d0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
