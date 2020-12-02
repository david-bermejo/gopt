#pragma once

#include "algorithms.hpp"
#include "constants.hpp"
#include "vector.hpp"
#include <cmath>

namespace gopt
{
	namespace internal_RKN7
	{
		template <typename T>
		static const T alpha[14] =
		{
			0,
			0.73470804841064383606103566874183e-1,
			0.11020620726159657540915535031127,
			0.16530931089239486311373302546691,
			0.5,
			0.26628826929126164263520369439706,
			0.63371173070873835736479630560294,
			0.75,
			0.5625,
			0.125,
			0.375,
			0.96652805160235700834451642119099,
			1,
			1
		};

		template <typename T>
		static const T beta[14][13] =
		{
			{0},
			{
				0.73470804841064383606103566874183e-1
			},
			{
				0.27551551815399143852288837577819e-1,
				0.82654655446197431556866512733456e-1
			},
			{
				0.41327327723098715778433256366728e-1,
				0,
				0.12398198316929614733529976910018
			},
			{
				0.89670558763795782658378325389875,
				0,
				-0.34585915268314336799411093327799e+1,
				0.30618859391934758533573260788811e+1
			},
			{
				0.57053369423965328229363646644660e-1,
				0,
				0,
				0.20664670695498245793175075857561,
				0.25881929123138564740892891767848e-2
			},
			{
				0.22130953402732738534298054285043e-1,
				0,
				0,
				0.36666901842159380713935315648106,
				0.32075560532767078702498038199180,
				-0.75843846443258975333835287154967e-1
			},
			{
				0.83333333333333333333333333333333e-1,
				0,
				0,
				0,
				0,
				0.38436436964131621037911008488971,
				0.28230229702535045628755658177696
			},
			{
				0.8349609375e-1,
				0,
				0,
				0,
				0,
				0.38306747479397408845870063561136,
				0.13548721270602591154129936438864,
				-0.3955078125e-1
			},
			{
				0.73420353223593964334705075445816e-1,
				0,
				0,
				0,
				0,
				0.98808964916022916024205710336420e-1,
				0.24153311327327749549842803451955,
				-0.48707561728395061728395061728395e-1,
				-0.24005486968449931412894375857339
			},
			{
				0.81378441127067064904041056404207e-2,
				0,
				0,
				0,
				0,
				0,
				-0.36266091174647134384031532058792,
				0.69726880597127928317272609847243e-1,
				0.37797780620763392161154341509711,
				0.28181838082900278742109519000315
			},
			{
				-0.14042538922482838913280031225476e+1,
				0,
				0,
				0,
				0,
				0.13555559029404957528304113342361e+2,
				-0,15021472824848050961721330969968e+1,
				0.14767543284167949686233606841588e+1,
				-0.21707681965133688432577373607995e+1,
				0.66149759502676558681039202833030e+1,
				0.11507526173569321530679222376434e+2
			},
			{
				-0.52708651815801315268176882187497e+1,
				0,
				0,
				0,
				0,
				-0.49965599553656833001045921529105e+2,
				-0.50302228928658231516135124812231e+1,
				0.44548269045298760506518238622704e+1,
				-0.86071533124033841312406742989148e+1,
				0.23840410046372287590078676456468e+2,
				0.41711581466028388124069667164840e+2,
				-0.13297747642437995408237095558512
			},
			{
				0.35099303056581883152660173681744e-1,
				0,
				0,
				0,
				0,
				0,
				0,
				0.25223475276631606400638853417712,
				0.11840033306876549234162515364336,
				0.20258133611250929893187899871888,
				0.26757025259420140796393329272621,
				0.16586384510629873791268098150965,
				-0.41749822704672884309167134456960e-1
			}
		};
		template <typename T>
		static const T gamma[14][13] =
		{
			{0},
			{
				0.26989795819968848329994970508715e-2
			},
			{
				0.30363520297464954371244341822304e-2,
				0.30363520297464954371244341822304e-2
			},
			{
				0.68317920669296147335299769100184e-2,
				0,
				0.68317920669296147335299769100184e-2
			},
			{
				-0.10263757731977888994310872824217e-2,
				0,
				0,
				0.12602637577319778889943108728242
			},
			{
				0.98909903843107417913313499064241e-2,
				0,
				0.20401758759111349514170518498571e-1,
				0.50265147713328703261825104735338e-2,
				0.13545726631277755415728360014730e-3
			},
			{
				0.36772464695317721429741572462201e-1,
				0,
				0,
				0.82132294778521785827721741407693e-1,
				0.30087165409098963036870918119641e-1,
				0.51803353935993790519824105531789e-1
			},
			{
				0.41233049088272873123221004021091e-1,
				0,
				0,
				0.11335100293061819105328798078376,
				0.56722148592237668841301677436715e-1,
				0.57456202064954525469376924736474e-1,
				0.12487597323916741512812413021961e-1
			},
			{
				0.4214630126953125e-1,
				0,
				0,
				0,
				-0.7808807373046875e-1,
				0.14104682102928772004397085536135,
				0.74603813736337279956029144638648e-1,
				-0.2150573730468750e-1
			},
			{
				0.55243877171925011431184270690444e-2,
				0,
				0,
				0,
				0,
				0.45913375893505158838018111807029e-2,
				0.12009956992268139808927955623138e-1,
				-0.24361818415637860082304526748971e-2,
				-0.11877000457247370827617741197988e-1
			},
			{
				0.12396099092300428073855581211091e-1,
				0,
				0,
				0,
				0,
				0,
				-0.23148568834881149606828637484335e-1,
				0.44057716338592294670599538200368e-2,
				0.24164236870396086181891828927171e-1,
				0.52494961238325405884021273526037e-1
			},
			{
				-0.12148292337172366838692371706654,
				0,
				0,
				-0.15948786809469047245658868763595e+1,
				0.77089844409590354601143580630038e-1,
				0,
				0,
				0.98844932135442618048624328441900e-1,
				-0.18517690177654009760124559303975,
				0.16665727117807342381867154630279e+1,
				0.52611925503652552568040599022802
			},
			{
				-0.49475846764102332689603709547782,
				0,
				0,
				-0.56513209641364305307648232070852e+1,
				0.42750028729043677987389324310306,
				0,
				0,
				0.30293416726956828108567954376506,
				-0.10280329379503342611614151091571e+1,
				0.54254171279669182157854764162220e+1,
				0.15340242607867031086671199895920e+1,
				-0.15763473585838266589893780962028e-1
			},
			{
				0.35173987586306713954725819037907e-1,
				0,
				0,
				0,
				0,
				0,
				0,
				0.63858784354258308506882892800822e-1,
				0.50866724905581448754291007148603e-1,
				0.17703179472766752427031494226269,
				0.16781715613041509463911067215393,
				0.45385629257942440722375392950401e-2,
				0.71298936997666580243712730101115e-3
			}
		};

		template <typename T, typename F, unsigned int N>
		void step(F& f, const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, Vector_t<Vector_t<T, N>, 14>& f_k, T t0, T h)
		{
			for (int k = 0; k < 14; k++)
			{
				Vector_t<T, N> x_k = x0 + dx0 * (alpha<T>[k] * h);
				Vector_t<T, N> dx_k = dx0;
				const T t_k = t0 + alpha<T>[k] * h;

				for (int i = 0; i < k; i++)
				{
					x_k += f_k[i] * (h * h * gamma<T>[k][i]);
					dx_k += f_k[i] * (h * beta<T>[k][i]);
				}

				f_k[k] = f(x_k, dx_k, t_k);
			}
		}

		template <typename T, unsigned int N>
		T calculate_error(const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, const Vector_t<Vector_t<T, N>, 14>& f_k, T h, T tol)
		{
			Vector_t<T, N> x_step = x0 + dx0 * h + (gamma<T>[13][0] * f_k[0] + gamma<T>[13][7] * f_k[7] + gamma<T>[13][8] * f_k[8] + gamma<T>[13][9] * f_k[9] + gamma<T>[13][10] * f_k[10] + gamma<T>[13][11] * f_k[11] + gamma<T>[13][12] * f_k[12]) * (h * h);
			Vector_t<T, N> rel_errors;

			for (int i = 0; i < N; i++)
				rel_errors[i] = std::abs((f_k[12][i] - f_k[13][i]) / x_step[i]);

			return max(rel_errors) * (gamma<T>[13][12] * h * h / tol);
		}

		template <typename T, unsigned int N>
		void update_step(Vector_t<T, N>& x_i, Vector_t<T, N>& dx_i, T& t_i, const Vector_t<Vector_t<T, N>, 14>& f_k, T h)
		{
			x_i += dx_i * h + (gamma<T>[13][0] * f_k[0] + gamma<T>[13][7] * f_k[7] + gamma<T>[13][8] * f_k[8] + gamma<T>[13][9] * f_k[9] + gamma<T>[13][10] * f_k[10] + gamma<T>[13][11] * f_k[11] + gamma<T>[13][12] * f_k[12]) * (h * h);
			dx_i += (beta<T>[13][0] * f_k[0] + beta<T>[13][7] * f_k[7] + beta<T>[13][8] * f_k[8] + beta<T>[13][9] * f_k[9] + beta<T>[13][10] * f_k[10] + beta<T>[13][11] * f_k[11] + beta<T>[13][12] * f_k[12]) * h;
			t_i += h;
		}
	}

	namespace internal_RKN5
	{
		template <typename T>
		static const T alpha[9] =
		{
			0,
			static_cast<T>(4) / 15,
			static_cast<T>(2) / 5,
			static_cast<T>(3) / 5,
			static_cast<T>(9) / 10,
			static_cast<T>(3) / 4,
			static_cast<T>(2) / 7,
			1,
			1
		};

		template <typename T>
		static const T beta[9][8] =
		{
			{0},
			{
				static_cast<T>(4) / 15
			},
			{
				static_cast<T>(1) / 10,
				static_cast<T>(3) / 10
			},
			{
				static_cast<T>(3) / 20,
				0,
				static_cast<T>(9) / 20
			},
			{
				static_cast<T>(9) / 40,
				0,
				0,
				static_cast<T>(27) / 40
			},
			{
				static_cast<T>(11) / 48,
				0,
				0,
				static_cast<T>(5) / 8,
				-static_cast<T>(5) / 48
			},
			{
				static_cast<T>(27112) / 194481,
				0,
				0,
				static_cast<T>(56450) / 64827,
				static_cast<T>(80000) / 194481,
				-static_cast<T>(24544) / 21609
			},
			{
				-static_cast<T>(26033) / 41796,
				0,
				0,
				-static_cast<T>(236575) / 38313,
				-static_cast<T>(14500) / 10449,
				static_cast<T>(275936) / 45279,
				static_cast<T>(228095) / 73788
			},
			{
				static_cast<T>(7) / 81,
				0,
				0,
				0,
				-static_cast<T>(250) / 3483,
				static_cast<T>(160) / 351,
				static_cast<T>(2401) / 5590,
				static_cast<T>(1) / 10
			}
		};

		template <typename T>
		static const T gamma[9][8] =
		{
			{0},
			{
				static_cast<T>(8) / 225
			},
			{
				static_cast<T>(1) / 25,
				static_cast<T>(1) / 25
			},
			{
				static_cast<T>(9) / 160,
				static_cast<T>(81) / 800,
				static_cast<T>(9) / 400
			},
			{
				static_cast<T>(81) / 640,
				0,
				static_cast<T>(729) / 3200,
				static_cast<T>(81) / 1600
			},
			{
				static_cast<T>(11283) / 88064,
				0,
				static_cast<T>(3159) / 88064,
				static_cast<T>(7275) / 44032,
				-static_cast<T>(33) / 688
			},
			{
				static_cast<T>(6250) / 194481,
				0,
				0,
				0,
				-static_cast<T>(3400) / 194481,
				static_cast<T>(1696) / 64827
			},
			{
				-static_cast<T>(6706) / 45279,
				0,
				0,
				0,
				static_cast<T>(1047925) / 1946997,
				-static_cast<T>(147544) / 196209,
				static_cast<T>(1615873) / 1874886
			},
			{
				static_cast<T>(31) / 360,
				0,
				0,
				0,
				0,
				static_cast<T>(64) / 585,
				static_cast<T>(2401) / 7800,
				-static_cast<T>(1) / 300
			}
		};

		template <typename T, typename F, unsigned int N>
		void step(F& f, const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, Vector_t<Vector_t<T, N>, 9>& f_k, T t0, T h)
		{
			for (int k = 0; k < 9; k++)
			{
				Vector_t<T, N> x_k = x0 + dx0 * (alpha<T>[k] * h);
				Vector_t<T, N> dx_k = dx0;
				T t_k = t0 + alpha<T>[k] * h;

				for (int i = 0; i < k; i++)
				{
					const T gamma_loc = gamma<T>[k][i];
					if (gamma_loc) x_k += f_k[i] * (h * h * gamma_loc);

					const T beta_loc = beta<T>[k][i];
					if (beta_loc) dx_k += f_k[i] * (h * beta_loc);
				}

				f_k[k] = f(x_k, dx_k, t_k);
			}
		}

		template <typename T, unsigned int N>
		T calculate_error(const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, const Vector_t<Vector_t<T, N>, 9>& f_k, T h, T tol)
		{
			Vector_t<T, N> x_step = x0 + dx0 * h + (gamma<T>[8][0] * f_k[0] + gamma<T>[8][5] * f_k[5] + gamma<T>[8][6] * f_k[6] + gamma<T>[8][7] * f_k[7]) * (h * h);
			Vector_t<T, N> rel_errors;

			for (int i = 0; i < N; i++)
				rel_errors[i] = std::abs((f_k[8][i] - f_k[7][i]) / x_step[i]);

			return max(rel_errors) * (h * h / (tol * 300));
		}

		template <typename T, unsigned int N>
		void update_step(Vector_t<T, N>& x_i, Vector_t<T, N>& dx_i, T& t_i, const Vector_t<Vector_t<T, N>, 9>& f_k, T h)
		{
			x_i += dx_i * h + (gamma<T>[8][0] * f_k[0] + gamma<T>[8][5] * f_k[5] + gamma<T>[8][6] * f_k[6] + gamma<T>[8][7] * f_k[7]) * (h * h);
			dx_i += (beta<T>[8][0] * f_k[0] + beta<T>[8][4] * f_k[4] + beta<T>[8][5] * f_k[5] + beta<T>[8][6] * f_k[6] + beta<T>[8][7] * f_k[7]) * h;
			t_i += h;
		}
	}

	template <typename T, typename F, unsigned int N>
	Vector_t<T, 2 * N> RKN7(F& f, const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, T t0, T tf, T tol)
	{
		T h = 0.001;
		T remaining = tf - t0;

		Vector_t<T, N> x_i = x0;
		Vector_t<T, N> dx_i = dx0;
		T t_i = t0;

		Vector_t<Vector_t<T, N>, 14> f_k;
		for (int i = 0; i < 14; i++)
			f_k[i] = Vector_t<T, N>(0);

		while (remaining > 0)
		{
			if (remaining < h)
			{
				h = remaining;
				internal_RKN7::step(f, x_i, dx_i, f_k, t_i, h);
				remaining = 0;

				internal_RKN7::update_step(x_i, dx_i, t_i, f_k, h);
				break;
			}

			internal_RKN7::step(f, x_i, dx_i, f_k, t_i, h);
			T error = internal_RKN7::calculate_error(x_i, dx_i, f_k, h, tol);
			T h_new = h;

			if (error < (static_cast<T>(1) / std::pow(2, 8)))
			{
				do
				{
					h_new *= 2;
					internal_RKN7::step(f, x_i, dx_i, f_k, t_i, h_new);
					error = internal_RKN7::calculate_error(x_i, dx_i, f_k, h_new, tol);
				} while (error < (static_cast<T>(1) / std::pow(2, 8)));

				if (error > 1)
				{
					h = h_new / 2;
					internal_RKN7::step(f, x_i, dx_i, f_k, t_i, h);
				}
				else
					h = h_new;
			}
			else if (error > 1)
			{
				do
				{
					h_new /= 2;
					internal_RKN7::step(f, x_i, dx_i, f_k, t_i, h_new);
					error = internal_RKN7::calculate_error(x_i, dx_i, f_k, h_new, tol);
				} while (error > 1);

				h = h_new;
			}

			internal_RKN7::update_step(x_i, dx_i, t_i, f_k, h);
			remaining -= h;
		}

		Vector_t<T, 2 * N> res;
		for (int i = 0; i < N; i++)
		{
			res[i] = x_i[i];
			res[N + i] = dx_i[i];
		}

		return res;
	}

	template <typename T, typename F, unsigned int N>
	Vector_t<T, 2 * N> RKN5(F& f, const Vector_t<T, N>& x0, const Vector_t<T, N>& dx0, T t0, T tf, T tol)
	{
		T h = 0.001;
		T remaining = tf - t0;

		Vector_t<T, N> x_i = x0;
		Vector_t<T, N> dx_i = dx0;
		T t_i = t0;

		Vector_t<Vector_t<T, N>, 9> f_k;
		for (int i = 0; i < 9; i++)
			f_k[i] = Vector_t<T, N>(0);

		while (remaining > 0)
		{
			if (remaining < h)
			{
				h = remaining;
				internal_RKN5::step(f, x_i, dx_i, f_k, t_i, h);
				remaining = 0;

				internal_RKN5::update_step(x_i, dx_i, t_i, f_k, h);
				break;
			}

			internal_RKN5::step(f, x_i, dx_i, f_k, t_i, h);
			T error = internal_RKN5::calculate_error(x_i, dx_i, f_k, h, tol);
			T h_new = h;

			if (error < (static_cast<T>(1) / std::pow(2, 6)))
			{
				do
				{
					h_new *= 2;
					internal_RKN5::step(f, x_i, dx_i, f_k, t_i, h_new);
					error = internal_RKN5::calculate_error(x_i, dx_i, f_k, h_new, tol);
				} while (error < (static_cast<T>(1) / std::pow(2, 6)));

				if (error > 1)
				{
					h = h_new / 2;
					internal_RKN5::step(f, x_i, dx_i, f_k, t_i, h);
				}
				else
					h = h_new;
			}
			else if (error > 1)
			{
				do
				{
					h_new /= 2;
					internal_RKN5::step(f, x_i, dx_i, f_k, t_i, h_new);
					error = internal_RKN5::calculate_error(x_i, dx_i, f_k, h_new, tol);
				} while (error > 1);

				h = h_new;
			}

			internal_RKN5::update_step(x_i, dx_i, t_i, f_k, h);
			remaining -= h;
		}

		Vector_t<T, 2 * N> res;
		for (int i = 0; i < N; i++)
		{
			res[i] = x_i[i];
			res[N + i] = dx_i[i];
		}

		return res;
	}

	namespace internal_DP45
	{
		template <typename T>
		static const T c[7] =
		{
			0,
			static_cast<T>(1) / 5,
			static_cast<T>(3) / 10,
			static_cast<T>(4) / 5,
			static_cast<T>(8) / 9,
			1,
			1
		};

		template <typename T>
		static const T a[7][6] =
		{
			{0},
			{
				static_cast<T>(1) / 5
			},
			{
				static_cast<T>(3) / 40,
				static_cast<T>(9) / 40
			},
			{
				static_cast<T>(44) / 45,
				-static_cast<T>(56) / 15,
				static_cast<T>(32) / 9
			},
			{
				static_cast<T>(19372) / 6561,
				-static_cast<T>(25360) / 2187,
				static_cast<T>(64448) / 6561,
				-static_cast<T>(212) / 729
			},
			{
				static_cast<T>(9017) / 3168,
				-static_cast<T>(355) / 33,
				static_cast<T>(46732) / 5247,
				static_cast<T>(49) / 176,
				-static_cast<T>(5103) / 18656
			},
			{
				static_cast<T>(35) / 384,
				0,
				static_cast<T>(500) / 1113,
				static_cast<T>(125) / 192,
				-static_cast<T>(2187) / 6784,
				static_cast<T>(11) / 84
			}
		};

		template <typename T>
		static const T b[2][7] =
		{
			{
				static_cast<T>(35) / 384,
				0,
				static_cast<T>(500) / 1113,
				static_cast<T>(125) / 192,
				-static_cast<T>(2187) / 6784,
				static_cast<T>(11) / 84,
				0
			},
			{
				static_cast<T>(5179) / 57600,
				0,
				static_cast<T>(7571) / 16695,
				static_cast<T>(393) / 640,
				-static_cast<T>(92097) / 339200,
				static_cast<T>(187) / 2100,
				static_cast<T>(1) / 40
			}
		};

		template <typename T>
		static const T e[7] =
		{
			b<T>[0][0] - b<T>[1][0],
			b<T>[0][1] - b<T>[1][1],
			b<T>[0][2] - b<T>[1][2],
			b<T>[0][3] - b<T>[1][3],
			b<T>[0][4] - b<T>[1][4],
			b<T>[0][5] - b<T>[1][5],
			b<T>[0][6] - b<T>[1][6]
		};

		template <typename T, typename F, unsigned int N>
		void step(F& f, const Vector_t<T, N>& x0, Vector_t<Vector_t<T, N>, 7>& f_k, T t0, T h)
		{
			for (int k = 0; k < 7; k++)
			{
				Vector_t<T, N> x_k = x0;
				T t_k = t0 + c<T>[k] * h;

				for (int i = 0; i < k; i++)
				{
					const T a_loc = a<T>[k][i];
					if (a_loc) x_k += f_k[i] * (h * a_loc);
				}

				f_k[k] = f(x_k, t_k);
			}
		}

		template <typename T, unsigned int N>
		T calculate_error(const Vector_t<T, N>& x0, Vector_t<T, N>& xnew, const Vector_t<Vector_t<T, N>, 7>& f_k, const T h, const T rtol, const T threshold)
		{
			Vector_t<T, N> fE = e<T>[0] * f_k[0];
			for (int i = 2; i < 7; i++)
				fE += e<T>[i] * f_k[i];
			
			xnew = x0 + (b<T>[0][0]*f_k[0] + b<T>[0][2]*f_k[2] + b<T>[0][3]*f_k[3] + b<T>[0][4]*f_k[4] + b<T>[0][5]*f_k[5]) * h;
			const Vector_t<T, N> xnew_abs = abs(xnew);
			const Vector_t<T, N> x_abs = abs(x0);

			Vector_t<T, N> tmp;
			for (int i = 0; i < N; i++)
				tmp[i] = fE[i] / std::max(std::max(x_abs[i], xnew_abs[i]), threshold);
			return max(abs(tmp)) * h;
		}

		template <typename T, unsigned int N>
		void update_step(Vector_t<T, N>& x_i, T& t_i, const Vector_t<Vector_t<T, N>, 7>& f_k, T h)
		{
			x_i += (b<T>[0][0]*f_k[0] + b<T>[0][2]*f_k[2] + b<T>[0][3]*f_k[3] + b<T>[0][4]*f_k[4] + b<T>[0][5]*f_k[5]) * h;
			t_i += h;
		}
	}

	template <typename T, typename F, unsigned int N>
	Vector_t<T,N> DP45(F& f, const Vector_t<T, N>& x0, const T t0, const T tf, const T rtol = static_cast<T>(1e-3), const T atol = static_cast<T>(1e-6), const T h0 = 0)
	{
		const T threshold = atol / rtol;
		T remaining = tf - t0;

		// Calculate initial h step if no one is provided.
		T h;
		if (!h0)
		{
			h = (tf - t0) / 10;
			const Vector_t<T, N> f0 = f(x0, t0);
			const Vector_t<T, N> x_abs = abs(x0);
			Vector_t<T, N> tmp;

			for (int i = 0; i < N; i++)
				tmp[i] = f0[i] / std::max(x_abs[i], threshold);

			const T rh = max(abs(tmp)) / (static_cast<T>(0.8) * std::pow(rtol, static_cast<T>(0.2)));

			if (h * rh > 1)
				h = 1 / rh;
		}
		else
			h = h0;

		Vector_t<T, N> xnew;

		// Initialize x_i and t_i variables.
		Vector_t<T, N> x_i = x0;
		T t_i = t0;

		// Declare the vector which holds the intermediate solutions of the function provided at each timestep.
		Vector_t<Vector_t<T, N>, 7> f_k;

		while (remaining > 0)
		{
			if (h * static_cast<T>(1.1) > remaining)
			{
				h = remaining;
				internal_DP45::step(f, x_i, f_k, t_i, h);
				remaining = 0;

				internal_DP45::update_step(x_i, t_i, f_k, h);
				break;
			}

			internal_DP45::step(f, x_i, f_k, t_i, h);
			const T error = internal_DP45::calculate_error(x_i, xnew, f_k, h, rtol, threshold);

			if (error > rtol)
				h *= std::max(static_cast<T>(0.1), static_cast<T>(0.8) * std::pow(rtol / error, static_cast<T>(0.2)));
			else
			{
				x_i = xnew;
				t_i += h;
				remaining -= h;

				const T tmp = static_cast<T>(1.25) * std::pow(error / rtol, static_cast<T>(0.2));
				h = (tmp > 0.2) ? h / tmp : h * 5;
			}
		}

		return x_i;
	}
}