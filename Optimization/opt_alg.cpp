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
		//std::cout << x0 << std::endl;
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
		double prev = X0.x(0);
		do
		{
			if (solution::f_calls > Nmax) {
				break;
				//throw "Nie znaleziono przedzialu po " + Nmax ;
			}
			
			prev = X0.x(0);
			X0 = X1;
			X1.x = x0 + alpha * d;
			X1.fit_fun(ff, ud1, ud2);
			alpha *= alpha;
			
		} while (X1.y <= X0.y);
		
		if (prev < X1.x) {
			p[0] = prev;
			p[1] = m2d(X1.x);
		}
		else
		{
			p[0] = m2d(X1.x);
			p[1] = prev;
		}

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
		//std::cout <<"DEPS:" << deps << std::endl;
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
		
		solution a0(a), b0(b);
		solution c0(b0.x - (static_cast<double>(tQ[k - 1]) / tQ[k]) * (b0.x - a0.x));
		solution d0(a0.x + b0.x - c0.x);
		Xopt.ud = b - a;
		for (int i = 1; i <= k - 3; i++)
		{
//std::cout << a0.x << b0.x;
//std::cout << c0.x << d0.x << std::endl;
			if (c0.fit_fun(ff, ud1, ud2) < d0.fit_fun(ff, ud1, ud2))
				b0 = d0;
			else
				a0 = c0;
			
			c0.x = b0.x - (static_cast<double>(tQ[k - i - 2]) / tQ[k - i - 1]) * (b0.x - a0.x);
			d0.x = a0.x + b0.x - c0.x;
//std::cout << a0.x << b0.x;
//std::cout << c0.x << d0.x << std::endl;
			Xopt.ud.add_row(m2d(b0.x - a0.x));
		}

		Xopt = c0;
		Xopt.flag = 0;
		//std::cout << Xopt.ud << std::endl << std::endl;
		delete[] tQ;
		
		return Xopt;
	}
	catch (string ex_info)
	{
		solution::clear_calls();
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		
		solution ai(a), bi(b), ci((a + b) * .5);
		ai.fit_fun(ff, ud1, ud2);
		bi.fit_fun(ff, ud1, ud2);
		ci.fit_fun(ff, ud1, ud2);
		long double di_1 = 0;
		int i = 0;
		solution di(0);
		Xopt.flag = 0;
		do
		{
			if (!i)
			{
				i++;
				Xopt.ud = b - a;
			}
			else
				Xopt.ud.add_row(m2d(bi.x - ai.x));
//std::cout << ai.x << bi.x << ci.x << std::endl;
//std::cout << ai.y << bi.y << ci.y << std::endl;

			matrix l = ai.y * (pow(bi.x, 2) - pow(ci.x, 2)) + bi.y * (pow(ci.x, 2) - pow(ai.x, 2)) 
				+ ci.y * (pow(ai.x, 2) - pow(bi.x, 2));
			matrix m = ai.y * (bi.x - ci.x) + bi.y * (ci.x - ai.x) + ci.y * (ai.x - bi.x);
//std::cout << m << l << std::endl;
			if (m <= 0)
			{
				Xopt.flag = -1;
				return Xopt;
			}

			di_1 = m2d(di.x);
			di = solution( 0.5 * m2d(l) / m2d(m));
//std::cout << di.x << std::endl << std::endl;
			if (abs(m2d(di.x) - di_1) < gamma)
				break;
			di.fit_fun(ff, ud1, ud1);

//std::cout << ai.x << " " << di.x << " " << bi.x << std::endl;
			if (ai.x < di.x && di.x < ci.x)
			{
				if (di.y < ci.y)
				{
					bi.x = ci.x;
					bi.y = ci.y;
					ci.x = di.x; 
					ci.y = di.y;

				}
				else
				{
					ai.x = di.x;
					ai.y = di.y;
				}

			}
			else if (ci.x < di.x && di.x < bi.x)
			{
				if (di.y < ci.y)
				{
					ai.x = ci.x;
					ai.y = ci.y;
					ci.x = di.x;
					ci.y = di.y;
				}
				else
				{
					bi.x = di.x;
					bi.y = di.y;
				}

			}
			else {
				Xopt.flag = -1;
				Xopt.x = di_1;
				Xopt.fit_fun(ff, ud1, ud2);
				return Xopt;
				//throw string("di poza zakresem\n");
			}
			if (solution::f_calls > Nmax) {
				
				throw string(string("Nie znaleziono przedzialu po " + Nmax) + " probach");
			}

		} while (m2d(bi.x) - m2d(ai.x) >= epsilon);
		
		Xopt.x = di.x;
		Xopt.y = di.y;
		
		Xopt.fit_fun(ff, ud1, ud2);
		// std::cout << Xopt.ud << std::endl << std::endl;
		
		return Xopt;
	}
	
	catch (string ex_info)
	{
		solution::clear_calls();
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution xs(x0);
		matrix ud(x0);

		xs.fit_fun(ff, ud1, ud2);
		Xopt.flag = 0;
		do
		{
			solution xB = xs;
			
			xs = HJ_trial(ff, xB, s, ud1, ud2); // check all directions from given base
			xs.fit_fun(ff, ud1, ud2);
			if (m2d(xs.y) < m2d(xB.y))
			{
				do
				{
					solution _xB = xB;
					xB = xs;
					xs = (xB.x * 2.0) - _xB.x; // calculate new point based on 2 last bases
					xs = HJ_trial(ff, xs, s, ud1, ud2); // check all directions from this new point
					if (solution::f_calls > Nmax)
					{
						Xopt.flag = -1;
						break;
						//throw string("Nie znaleziono przedzialu po Nmax probach (f(x)<f(xB))");
					}
				} while (m2d(xs.y) < m2d(xB.y)); // end this loop when the algorithm couldn't find better base with given step
				xs = xB;
				
			}
			else
			{
				// if the search function is unable to find any better point in proximity of the base, make step smaller
				s = alpha * s;
			}
			if (solution::f_calls > Nmax)
			{
				if(!Xopt.flag)
					Xopt.flag = -2;
				break;
				//throw string("Nie znaleziono przedzialu po Nmax probach (L)");
			}
			ud.add_col(xs.x);
			// if step size is smaller than targeted precision, return the last calculated point
		} while (s >= epsilon);
		
		Xopt = xs;
		Xopt.ud = ud;
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
		const int DIM = 2;
		
		matrix dj(DIM, DIM);
		
		for(int k = 0; k < DIM; k++)
			for (int l = 0; l < DIM; l++)
				dj(k, l) = (l == k)? 1 : 0;
			
		for (int j = 0; j < DIM; j++)
		{
			// calculate new point, if new value is closer to minimum, change point to be searched 
			solution xj = XB.x + (s * dj[j]); 
			xj.fit_fun(ff, ud1, ud2);
			if ( xj.y < XB.y ) 
			{
				XB = xj;
			}
			else
			{
				// if not, check if the mirror point is better
				xj = XB.x - (s * dj[j]);
				xj.fit_fun(ff, ud1, ud2);
				if (xj.y < XB.y)
				{
					XB = xj;
				}
			}
		}

		
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
		
		const int DIM = 2;
		//matrix: row column
		matrix d(DIM, DIM); // DIM x DIM square matrix
		matrix ud(x0);

		for (int w = 0; w < DIM; w++)
			for (int k = 0; k < DIM; k++)
				d(w, k) = (w == k) ? 1 : 0;

		matrix l(DIM, 1, 0.0); // vertical vector
		matrix p(DIM, 1, 0.0); // vertical vector non-essential
		matrix s(s0); // vertical vector

		solution xB(x0); // x -> vertical vector; y -> scalar
		xB.fit_fun(ff, ud1, ud2);
		solution Xopt(xB);
		Xopt.flag = 0;
		int max_s;
		do
		{
			for (int j = 0; j < DIM; j++)
			{
				solution _x(xB.x + s(j) * d[j]);
				if (_x.fit_fun(ff, ud1, ud2) <  xB.y) // better solution, change base, extend searching reach
				{
					xB = _x;
					
					l(j) = l(j) + s(j);
					s(j) = s(j) * alpha;
					
				}
				else // report the failure, change step (reverse it and make it bigger)
				{
					
					s(j) = -s(j) * beta;
					p(j) = p(j) + 1;
					
				}
			}
			// after calculating new base point, replace the old answer with it
			Xopt = xB;

			// check, if there was at least one success and one failure on each orthonormal direction
			bool zero = false;
			for (int j = 0; j < DIM; j++)
			{
				if (p(j) == 0 || abs(l(j)) < epsilon ) { // compare double with given precision
					zero = true; // even if only one direction have not been examied, do not make any changes
					break;
				}
			}

			if (!zero)
			{
				// change direction { square matrix d(DIM, DIM) }
				matrix _D(d); // square matrix initialised with d matrix values
				matrix _lQ(DIM, DIM); // empty square matrix (DIM, DIM)

				// create triangle matrix based on lambda values
				for (int i = 0; i < DIM; i++) // row
					for (int j = 0; j < DIM; j++) // column
						_lQ(i, j) = (i >= j) ? l(i) : 0.0;

				// calulate Q square matrix (DIM, DIM)
				_lQ = _D * _lQ;

				// create and set new directional values, starting from first column (first direction vector)
				matrix v(DIM,DIM);
				v.set_col(_lQ[0]/(norm(_lQ[0])), 0);
				
				for (int _j = 1; _j < DIM; _j++) // v/Q matrix column
				{
					matrix sigma(DIM,1);
					matrix t_lQ(trans(_lQ[_j]));

					for (int k = 0; k < _j; k++) // sigm column, ie. current dimension which have it's orientation calculated, all vectors must be orthogonal
					{
						sigma.set_col(
							sigma[0] + (t_lQ * d[k]) * d[k],
							0);
						
					}
					// pk represents vertical vector, which must be normalised 
					matrix pk = _lQ[_j] - sigma[0];
					v.set_col(pk/norm(pk), _j);
					
				}
				
				// override previous direction(orientation) matrix with the new one
				d = v;
				
				// reverse lambda, failure and step matrixes to initial values
				l = matrix(DIM, 1, 0.0);
				p = matrix(DIM, 1, 0.0);
				s = s0;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
				//throw string("Nie znaleziono przedzialu po Nmax probach (f_calls > Nmax)");
			}
			ud.add_col(Xopt.x);
			// find the largest step and compare it to the targeted precision
			max_s = 0;
			for (int j = 1; j < DIM; j++)
			{
				if (abs(s(max_s)) < abs(s(j)))
				{
					max_s = j;
				}
			}
			
		} while (abs(s(max_s)) >= epsilon);
		Xopt.ud = ud;
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
		solution xi(x0), x_i;
		xi.flag = 0;
		matrix init_v_S(2, 1);
		init_v_S(0) = c;
		init_v_S(1) = ud2(0);
		double nm = 0;
		//double tmp1 = m2d(ff3T(xi.x)), tmp2;
		do
		{
			x_i = xi;
			//tmp2 = tmp1;
			xi = sym_NM(ff,xi.x,ud1(0),ud1(1),ud1(2),ud1(3), ud1(4), ud1(5), Nmax, init_v_S);
			//tmp1 = m2d(ff3T(xi.x));
//std::cout << "PEN:\nx:" << x_i.x(0) << " " << x_i.x(1) << " y: " << xi.y << "\nx:" << xi.x(0) << " " << xi.x(1) << " y: " << xi.y << "\n";
			init_v_S(0) = init_v_S(0) * dc;
			if (solution::f_calls > Nmax)
			{
				xi.flag = -2;
				break;
			}
			if (dc < 1.0)
			{
				double sum = 0.0;
				sum += 1 / m2d(g1(xi.x, NAN));
				sum += 1 / m2d(g2(xi.x, NAN));
				sum += 1 / m2d(g3(xi.x, init_v_S(1)));
				if (c * fabs(sum) < epsilon)
					break;
			}
			
			nm = norm(xi.x - x_i.x);

		} while (nm >= epsilon);
		
		return xi;
	}
	catch (string ex_info)
	{
		throw ("matrix pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		
		solution Xopt;
		int g_min = 0;
		Xopt.flag = 0;
		// set function's dimension
		const int DIM = 2;
		// set initial matrixes
		matrix p(DIM, 1 + DIM); // DIM (points' dimension) x DIM+1 (points' amount) square point matrix
		matrix e = ident_mat(DIM); // DIM x DIM square identity matrix
		
		p.set_col(x0, 0); // p0 = x0
		for (int i = 1; i <= DIM; i++)
			p.set_col(p[0] + e[i-1] * s, i); //pi = p0 + ei*s

		solution::f_calls += 1 + DIM;
		matrix p_f(DIM+1, 1);
		for (int i = 0; i <= DIM; i++)
			p_f(i) = m2d(ff(p[i], ud1, NAN)); // returns matrix 1x1
		
		double max_norm;
		do {
			max_norm = 0.0;
			int p_max = 0, p_min = 0;
			for (int i = 1; i <= DIM; i++) {
				if (p_f(p_max) < p_f(i)) p_max = i;
				if (p_f(p_min) > p_f(i)) p_min = i;
			}
			if (p_max == p_min)
				p_max = (p_max+1)%(DIM+1);

			matrix p_s(DIM,1); // _p
			for (int i = 0; i <= DIM; i++)
			{
				if (i == p_max) continue;
				p_s.set_col(p_s[0] + p[i], 0); // _p = E(i!=max) pi
			}

			p_s.set_col(p_s[0] / DIM, 0);  // _p /= n
			matrix p_odb = p_s[0] + (p_s[0] - p[p_max]) * alpha; // p_odb = _p + a(_p - p_max)
			
			solution::f_calls++;
			double p_odb_f = m2d(ff(p_odb, ud1, NAN)); // returns matrix 1x1
			
			if (m2d(p_odb_f) < p_f(p_max))
			{
				matrix p_e = p_s + (p_odb[0] - p_s[0]) * gamma;

				solution::f_calls++;
				double p_e_f = m2d(ff(p_e, ud1, NAN));

				if (ff(p_e, ud1, NAN) < p_odb_f)
				{
					p.set_col(p_e[0], p_max);
					p_f(p_max) = p_e_f;
				}
				else
				{
					p.set_col(p_odb[0], p_max);
					p_f(p_max) = p_odb_f;
				}
			}
			else
			{
				if (p_f(p_min) <= p_odb_f && p_odb_f < p_f(p_max))
				{
					p.set_col(p_odb[0], p_max);
					p_f(p_max) = p_odb(0);
				}
				else
				{
					matrix p_z = p_s[0] + (p[p_max] - p_s[0])*beta;

					solution::f_calls++;
					double p_z_f = m2d(ff(p_z, ud1, NAN));

					if (p_z_f >= p_f(p_max))
					{
						for (int i = 0; i <= DIM; i++)
						{
							if (i == p_min) continue;
							p.set_col((p[i] + p[p_min]) * delta, i);
							solution::f_calls++;
							p_f(i) = m2d(ff(p[i], ud1, NAN));
						}
					}
					else
					{
						p.set_col(p_z[0], p_max);
						p_f(p_max) = p_z_f;
					}
				}
			}
			g_min = p_min;
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
				//throw("Nie znaleziono przedzialu po Nmax probach (f_calls > Nmax)\n");
			}
			for (int i = 0; i <= DIM; i++)
			{
				if (i == p_min) continue;
				double i_norm = norm(p[p_min] - p[i]);
				if (i_norm > max_norm)
					max_norm = i_norm;
			}
			//std::cout << max_norm << std::endl;
		} while (max_norm >= epsilon);
		Xopt.x = p[g_min];
		Xopt.y = p_f(g_min);
//std::cout << "NM: " << Xopt << std::endl;
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
		solution xi(x0), x_i;
		xi.flag = 0;
		bool h0_c = false;
		if (h0 == 0.0)
			h0_c = true;
		int i = 0;
		bool save_prev_x = !isnan(ud1(0));
		
		do
		{
			matrix di = -gf(xi.x, NAN, NAN);
			++solution::g_calls;
			if (h0_c) {
				double* range = expansion(ff, 0, 10.0, 2.0, Nmax, xi.x, di);
				h0 = m2d(golden(ff, 0, range[1], epsilon, Nmax, xi.x, di).x);
				delete[] range;
			}

			x_i = xi;
			if (save_prev_x)
			{
				if (!i) {
					xi.ud = matrix(x_i.x);
					i++;
				}
				else {
					xi.ud.add_col(x_i.x);
				}
			}
			xi.x = xi.x + di * h0;
			
			xi.fit_fun(ff,NAN,NAN);

			if (solution::f_calls > Nmax) {
				xi.flag = -2;
				break;
			}
			
		} while (norm(xi.x - x_i.x) >= epsilon);

		return xi;
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
		solution xi(x0), x_i;
		xi.flag = 0;
		bool h0_c = false;
		if (h0 == 0.0)
			h0_c = true;
		bool save_prev_x = !isnan(ud1(0));
		int i = 0;
		int k = 0;
		matrix di;
		double x_i_g_pow_norm =0.0;
		do
		{
			// Debug output
			if (solution::f_calls % 10000 == 0 && solution::f_calls != 0)
				std::cout << "|" << solution::f_calls / 10000 << "|";

			xi.grad(gf, ud1, ud2);
			double xi_g_pow_norm = pow(norm(xi.g), 2);
			
			if (i != 0) {

				double beta = xi_g_pow_norm / x_i_g_pow_norm;
				if (x_i.y < xi.y)
				{
					//std::cout << "i: " << i << " " << x0 << std::endl;
					di = -x_i.g;
				}
				else 
				di = -xi.g + di * beta;
			}
			else {
				di = -xi.g;
				
			}
			i++;
			
			x_i_g_pow_norm = xi_g_pow_norm;
			
			if (h0_c) {
				double* range = expansion(ff, 0, 10.0, 2.0, Nmax, xi.x, di);
				h0 = m2d(golden(ff, max(range[0],0.0), range[1], epsilon, Nmax, xi.x, di).x);
				delete[] range;
			}
			
			x_i.x = xi.x;
			x_i.y = xi.y;
			x_i.g = xi.g;
			if (save_prev_x)
			{
				if (!k) {
					xi.ud = matrix(x_i.x);
					k++;
				}
				else {
					xi.ud.add_col(x_i.x);
				}
			}

			xi.x = xi.x + di * h0;
			xi.fit_fun(ff, ud1, ud2);
			

			if (solution::f_calls > Nmax) {
				xi.flag = -2;
				break;
			}

		} while (norm(xi.x - x_i.x) >= epsilon);

		return xi;
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
		solution xi(x0), x_i;
		xi.flag = 0;
		bool h0_c = false;
		if (h0 == 0.0) 
			h0_c = true;
		bool save_prev_x = !isnan(ud1(0));
		int i = 0;
		do
		{
			xi.hess(Hf, NAN, NAN);
			xi.grad(gf, NAN, NAN);
			
			matrix di = -inv(xi.H)*xi.g;

			if (h0_c) {
				double* range = expansion(ff, 0, 10.0, 2.0, Nmax, xi.x, di);
				h0 = m2d(golden(ff, max(range[0], 0.0), range[1], epsilon, Nmax, xi.x, di).x);
				delete[] range;
			}

			if (solution::f_calls > Nmax)
			{
				xi.flag = -2;
				break;
			}

			x_i = xi;
			if (save_prev_x)
			{
				if (!i) {
					xi.ud = matrix(xi.x);
					
					i++;
				}
				else {
					xi.ud.add_col(xi.x);
				}
			}

			xi.x = xi.x + di * h0;

			xi.fit_fun(ff, NAN, NAN);

			solution::f_calls++;
			if (solution::f_calls > Nmax) {
				xi.flag = -2;
				break;
			}

		} while (norm(xi.x - x_i.x) >= epsilon);

		return xi;
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
		Xopt.flag = 0;
		//Tu wpisz kod funkcji
		static const double alpha = (sqrt(5.0) - 1.0) * .5;
		double ai = a, bi = b;
		double ay = m2d(ff(a, ud1, ud2)), by = m2d(ff(b, ud1, ud2));
		solution::f_calls+= 2;
		double ci  = (bi - alpha * (bi - ai)), di = (ai + alpha * (bi - ai));
		double cy = m2d(ff(ci, ud1, ud2)), dy = m2d(ff(di, ud1, ud2));
		solution::f_calls += 2;
		do
		{
			if (cy < dy)
			{
				bi = di;
				by = dy;
				di = ci;
				dy = cy;
				ci = bi - alpha * (bi - ai);
				cy = m2d(ff(ci, ud1, ud2));
				solution::f_calls++;
			}
			else
			{
				ai = ci;
				ay = cy;
				ci = di;
				cy = dy;
				di = ai + alpha * (bi - ai);
				dy = m2d(ff(di, ud1, ud2));
				solution::f_calls++;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -2;
				break;
			}
		} while (bi-ai >= epsilon);
		Xopt.x = matrix((ai + bi) / 2.0);
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
		int DIM = get_len(x0);
		matrix e = ident_mat(DIM);
		matrix d = matrix(e);
		matrix xi = x0;
		do
		{
			matrix p = matrix(DIM,DIM+1);
			p.set_col(xi,0);
			for(int j = 1; j <= DIM; j++)
			{
				//wyznacz h
				//std::cout << " range: \n";
				double* range = expansion(ff, -100.0, 100.0, 2.0, Nmax, p[j - 1], d[j-1]);
				double h0 = m2d(golden(ff, range[0], range[1], epsilon, Nmax, p[j - 1], d[j - 1]).x);
				p.set_col(p[j - 1] + d[j - 1] * h0,j);
				delete[] range;
			}
			if (norm(p[DIM] - xi) < epsilon)
			{
				Xopt.x = xi;
				Xopt.flag = 0;
				break;
			}
			for (int j = 0; j < DIM-1; j++)
			{
				d.set_col(d[j+1],j);
			}
			d.set_col(p[DIM] - p[0], DIM - 1);

			double* range = expansion(ff, 0, 10.0, 2.0, Nmax, p[DIM], d[DIM - 1]);
			double h0 = m2d(golden(ff, range[0], range[1], epsilon, Nmax, p[DIM], d[DIM - 1]).x);
			xi = p[DIM] + h0 * d[DIM - 1];
			delete[] range;
			if (solution::f_calls > Nmax)
			{
				Xopt.x = xi;
				Xopt.flag = -2;
				break;
			}
		} while(true);
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
