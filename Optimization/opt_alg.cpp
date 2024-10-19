#include"opt_alg.h"

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
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji
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
				p[1] = m2d(X0.x)-d;
				
				return p;
			}
		}
		double prev;
		do
		{
			
			if (solution::f_calls > Nmax) {
				solution::clear_calls();
				throw string(string("Nie znaleziono przedzialu po " +  Nmax) + " probach");
			}
			prev = m2d(X0.x);
			X0 = X1;
			X1.x = x0 + alpha * d;
			X1.fit_fun(ff, ud1, ud2);
			alpha *= alpha;
			
		} while (X1.y <= X0.y);
		if (d > 0)
		{
			p[0] = prev;
			p[1] = m2d(X1.x);
		}
		else
		{
			p[1] = m2d(X1.x);
			p[0] = prev;
		}
		solution::clear_calls();
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
		//Tu wpisz kod funkcji
		 
		//znajdz najmniejsza k, tak aby k-ty element ciagu fibonaciego byl wiekszy od dlugosci przedzialu przez epsylon
		double deps = (b - a) / epsilon;
		std::cout <<"DEPS:" << deps << std::endl;
		int k = 2;
		int Mk = 10;
		const int ak = 5;
		int Pk = 10;
		long* tQ = new long[Mk];
		tQ[0] = 0;
		tQ[1] = 1;
		tQ[2] = 1;
		tQ[3] = 2;
		tQ[4] = 3;
		tQ[5] = 5;
		tQ[6] = 8;
		tQ[7] = 13;
		tQ[8] = 21;
		tQ[9] = 34;
		while (tQ[k] <= deps)
		{
			k++;
			
			if (k != Mk)
			{
				Mk += ak;
				long* tmpQ = new long[Mk];
				std::copy(tQ, tQ + Mk - ak, tmpQ);
				delete[] tQ;
				tQ = tmpQ;
			}

			if (k == Pk) {
				
				tQ[k] = tQ[k - 1] + tQ[k - 2];
				Pk++;
			}

		}
		
		double a0 = a, b0 = b;
		solution c0(b0 - static_cast<double>(tQ[k - 1]) / tQ[k] * (b0 - a0));
		solution d0(static_cast<double>(tQ[k - 1]) / tQ[k] * (b0 - a0));
		for (int i = 0; i <= k - 3; i++)
		{

			if (c0.fit_fun(ff, ud1, ud2) < d0.fit_fun(ff, ud1, ud2))
				b0 = m2d(d0.x);
			else
				a0 = m2d(c0.x);

			c0.x = b0 - static_cast<double>(tQ[k - i - 2]) / tQ[k - i - 1] * (b0 - a0);
			d0.x = a0 + b0 - c0.x;
		}

		Xopt = c0;
		Xopt.flag = 0;
		return Xopt;
		delete[] tQ;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		int i = 0;
		solution ai(a), bi(b), ci((a + b) * .5);
		ai.fit_fun(ff, ud1, ud2);
		bi.fit_fun(ff, ud1, ud2);
		ci.fit_fun(ff, ud1, ud2);
		double di_1 = 0;
		solution di(0);
		do
		{
			matrix l = ai.y * (pow(bi.x, 2) - pow(ci.x,2)) + bi.y * (pow(ci.x, 2) - pow(ai.x, 2)) 
				+ ci.y * (pow(ai.x, 2) - pow(bi.x, 2));
			matrix m = ai.y * (bi.x - ci.x) + bi.y * (bi.x - ci.x) + ci.y * (ai.x - bi.x);

			if (m <= 0)
				throw string("Brak rozwi¹zania\n");

			di_1 = m2d(di.x);
			di = solution( 0.5 * m2d(l) / m2d(m));
			di.fit_fun(ff, ud1, ud1);

			if (abs(m2d(di.x) - di_1) < gamma)
				break;
			if (ai.x < di.x && di.x < ci.x)
			{
				if (di.y < ci.y)
				{
					a = m2d(ai.x);
					ci = di;
					bi = ci;
				}
				else
				{
					ai = di;
					b = m2d(bi.x);
				}
			}
			else if (ci.x < di.x && di.x < bi.x)
			{
				if (di.y < ci.y)
				{
					a = m2d(ci.x);
					ai = ci;
					ci = di;
				}
				else
				{
					a = m2d(ai.x);
					bi = di;
				}
			}
			else
				throw string("d(i) poza zakresem\n");
			i++;
			if (solution::f_calls > Nmax) {
				throw string(string("Nie znaleziono przedzialu po " + Nmax) + " probach");
			}

		} while (!(b-a < epsilon));
		Xopt = di;
		Xopt.flag = 0;
		
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
