#include <iostream>
#include "math/gopt.hpp"
using namespace gopt;

#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;

struct State
{
	union {
		struct {
			Vec3 x;
			Vec3 v;
			Vec3 w;
			Quat q;
		};

		Vector_t<double, 13> value;
	};

	State() : value(0) {}
};

struct Params
{
	double m0;
	double m_dot;
	double Ixx;
	double Iyy;
	double Izz;
	double Ixx_dot;
	double Iyy_dot;
	double Izz_dot;
	double d_com_ct;
	Vec2 tvc;
	double tb;
	double d_com_cp;
	double Cd0;
	double Cda;
	double Cla;
	double rho;
	double Sref;
	double T;
};

Vector_t<double, 13> dynamics(const Vector_t<double, 13>& x, const double& t, const Params& p)
{
	State& s = *(State*)(&x[0]);

	const double T_val = (t < p.tb) ? p.T : 0;

	const double sx = std::sin(p.tvc[0]);
	const double cx = std::cos(p.tvc[0]);
	const double sy = std::sin(p.tvc[1]);
	const double cy = std::cos(p.tvc[1]);
	const Vec3 T = Vec3(sy, -sx * cy, cx * cy) * T_val;

	const double m = p.m0 + p.m_dot * std::min(t, p.tb);
	const double Ixx = p.Ixx + p.Ixx_dot * std::min(t, p.tb);
	const double Iyy = p.Iyy + p.Iyy_dot * std::min(t, p.tb);
	const double Izz = p.Izz + p.Izz_dot * std::min(t, p.tb);

	// Aerodynamic forces
	const double v_mag = s.v.length();
	Vec3 F_ae;
	double alpha = 0;

	if (v_mag > 1e-6)
	{
		const Vec3 nz(0, 0, 1);
		const Vec3 vu = rotate(s.q, s.v / v_mag);
		const Vec3 nD = -vu;
		Vec3 nL = nz - vu[2]*vu;
		nL = (nL.length() > 1e-6) ? normalize(nL) : Vec3(0);
		alpha = std::abs(std::acos(std::clamp(vu[2], -1.0, 1.0)));

		const double rho = 1.225;
		const Vec3 L = 0.5 * rho * v_mag * v_mag * p.Cla * alpha * p.Sref * nL;
		const Vec3 D = 0.5 * rho * v_mag * v_mag * (p.Cd0 + p.Cda * alpha) * p.Sref * nD;
		F_ae = L + D;
	}
	else
		F_ae = Vec3(0);

	const Vec3 forces = rotate(conjugate(s.q), T + F_ae) + Vec3(0, 0, -m * 9.81);
	const Vec3 moments = cross(Vec3(0, 0, -p.d_com_ct), T) + cross(Vec3(0, 0, p.d_com_cp), F_ae);

	const Mat3 I(Ixx, 0, 0, 0, Iyy, 0, 0, 0, Izz);
	const Mat3 I_inv(1/Ixx, 0, 0, 0, 1/Iyy, 0, 0, 0, 1/Izz);
	const Mat3 I_dot = (t < p.tb) ? Mat3(p.Ixx_dot, 0, 0, 0, p.Iyy_dot, 0, 0, 0, p.Izz_dot) : Mat3(0);

	State res;
	res.x = s.v;
	res.v = forces / m;
	res.w = I_inv * (moments - I_dot * s.w - cross(s.w, I * s.w));
	res.q = -(Quat(0.0, s.w[0], s.w[1], s.w[2]) * normalize(s.q)) / 2;

	return res.value;
}

Vec2 project(Quat q, Vec3 axis)
{
	const Vec3 v = rotate(q, axis);
	return { std::atan2(-v[1], v[2]), std::atan2(v[0], v[2]) };
}

struct Result
{
	std::vector<double> x;
	std::vector<Vector_t<double, 19>> y;
};

double objective(const Vector_t<double, 6>& K, const double phi, const double theta, const double psi, const double tvc_bound, const bool print = false, Result* result = nullptr)
{
	double prev_t = 0;
	const double dt = 0.01;

	// This is right
	State state;

	Quat qx = Quat(phi, Vec3(1, 0, 0));
	Quat qy = Quat(theta, Vec3(0, 1, 0));
	Quat qz = Quat(psi, Vec3(0, 0, 1));

	state.q = qx * qy * qz;

	Params params;
	params.m0 = 0.85; // kg
	params.tb = 5; // s
	params.m_dot = -0.2 * params.m0 / params.tb; // kg/s

	params.Ixx = 0.05;
	params.Iyy = 0.05;
	params.Izz = 0.01;
	params.Ixx_dot = -0.2 * params.Ixx / params.tb;
	params.Iyy_dot = -0.2 * params.Iyy / params.tb;
	params.Izz_dot = -0.2 * params.Izz / params.tb;

	params.d_com_ct = 0.4; // m
	params.tvc = 0; // rad
	params.d_com_cp = 0.2; // m
	params.Cd0 = 0.15;
	params.Cda = 0.25;
	params.Cla = 0.1;
	params.rho = 1.225;
	params.Sref = pi<> * 0.045 * 0.045;
	params.T = 11; // N

	const Mat2 Kp(K[0], 0, 0, K[1]);
	const Mat2 Ki(K[2], 0, 0, K[3]);
	const Mat2 Kd(K[4], 0, 0, K[5]);

	Vec2 error = 0;
	Vec2 prev_error = 0;
	Vec2 prev2_error = 0;
	Vec2 int_accum = 0;

	double D = 0;
	double t = 0;

	while (true)
	{
		if (prev_t < params.tb)
		{
			// PID
			const Vec2 set_point = 0;
			error = set_point - project(state.q, Vec3(0, 0, 1));

			// Proportional
			Vec2 control = Kp * error;

			// Integral
			int_accum += error + prev_error;
			control += Ki * int_accum * (dt/2);

			// Derivative (2nd order aproximation)
			control += Kd * (3 * error - 4 * prev_error + prev2_error) / (2 * dt);

			// Transform forces into tvc angles
			params.tvc[0] = std::clamp(control[0], -tvc_bound, tvc_bound);
			params.tvc[1] = std::clamp(control[1], -tvc_bound, tvc_bound);

			const double cx = std::cos(params.tvc[0]);
			const double cy = std::cos(params.tvc[1]);
			const double tz = cx * cy;

			if (std::acos(tz) > tvc_bound)
			{
				const double sx = std::sin(params.tvc[0]);
				const double sy = std::sin(params.tvc[1]);

				const double tx = -sy;
				const double ty = -sx * cy;

				const double a = std::pow(std::tan(std::atan2(ty, tx)), 2);
				const double b = std::pow(std::cos(tvc_bound), 2);

				params.tvc[1] = copysign(std::asin(std::sqrt((1 - b) / (1 + a))), params.tvc[1]);
				params.tvc[0] = copysign(std::acos(std::sqrt(b) / std::cos(params.tvc[1])), params.tvc[0]);
			}

			// Update previous variables
			prev2_error = prev_error;
			prev_error = error;

			D += dt / 2 * (error.length() + prev_error.length());
		}
		else
			params.tvc = 0;

		// System propagation
		t = prev_t + dt;
		const auto fcn = [&](const Vector_t<double, 13> x, double t) { return dynamics(x, t, params); };
		const auto res = DP45(fcn, state.value, prev_t, t, 1e-9, 1e-9);

		state.value = res;
		prev_t = t;

		if (result)
		{
			result->x.push_back(t);

			Vector_t<double, 19> v;
			for (int i = 0; i < 13; i++)
				v[i] = res[i];

			for (int i = 0; i < 2; i++)
				v[13 + i] = error[i];

			v[15] = params.tvc[0];
			v[16] = params.tvc[1];

			// Dynamic pressure
			v[17] = 0.5 * 1.225 * state.v.magnitude();

			// Angle of attack
			v[18] = state.v.magnitude() ? std::abs(std::acos(std::clamp(normalize(rotate(state.q, state.v))[2], -1.0, 1.0))) : 0;

			result->y.push_back(v);
		}

		if ((t > params.tb) && state.x[2] <= 0 || state.v[2] <= 0)
			break;
	}

	D -= t;
	if (print)
		printf("objective: %g, -> [Kp: %g %g, Ki: %g %g, Kd: %g %g], t_final: %g\n", D + t, K[0], K[1], K[2], K[3], K[4], K[5], t);

	return D;
}

int main()
{
	const double Kp = 10;
	const double Ki = 1;
	const double Kd = 1;

	const Vector_t<double, 6> ub(Kp,Kp, Ki,Ki, Kd,Kd);
	const Vector_t<double, 6> lb = -ub;

	// Gravity vector
	const double alpha = radians(245);
	const double delta = radians(8.0);
	const double R = std::tan(delta) / (1 + std::tan(delta) * std::tan(delta));
	Vec3 grav = Vec3(R * std::cos(alpha), R * std::sin(alpha), std::sqrt(1 - R * R));

	const double phi = std::atan2(grav[1], grav[2]);
	const double theta = std::atan2(-grav[0], std::sqrt(grav[1] * grav[1] + grav[2] * grav[2]));

	const double tvc_bound = radians(10.0);

	const auto fcn = [&](const Vector_t<double, 6>& x)
	{
		return objective(x, phi, theta, 0.0, tvc_bound);
	};

	const auto best = particleswarm(fcn, lb, ub, 100, 10);
	std::cout << "Best: " << best << ", objective value: " << std::endl;

	Result result;
	objective(best, phi, theta, 0, tvc_bound, true, &result);

	std::vector<double> x(result.y.size());
	std::vector<double> y(result.y.size());
	std::vector<double> z(result.y.size());

	std::vector<double> vel(result.y.size());

	std::vector<double> err_x(result.y.size());
	std::vector<double> err_y(result.y.size());

	std::vector<double> theta_x(result.y.size());
	std::vector<double> theta_y(result.y.size());

	std::vector<double> dyn_pressure(result.y.size());
	std::vector<double> AoA(result.y.size());

	for (int i = 0; i < result.x.size(); i++)
	{
		x[i] = result.y[i][0];
		y[i] = result.y[i][1];
		z[i] = result.y[i][2];

		vel[i] = std::sqrt(std::pow(result.y[i][3], 2) + std::pow(result.y[i][4], 2) + std::pow(result.y[i][5], 2));

		err_x[i] = degrees(result.y[i][13]);
		err_y[i] = degrees(result.y[i][14]);

		theta_x[i] = degrees(result.y[i][15]);
		theta_y[i] = degrees(result.y[i][16]);

		dyn_pressure[i] = result.y[i][17];
		AoA[i] = degrees(result.y[i][18]);
	}

	// Plot
	plt::figure();
	plt::named_plot("x", result.x, x, "-");
	plt::named_plot("y", result.x, y, "-");
	plt::named_plot("z", result.x, z, "-");
	plt::title("Position");
	plt::xlabel("t (s)");
	plt::ylabel("position (m)");
	plt::legend();
	plt::xlim(0.0, result.x.back());

	plt::figure();
	plt::named_plot("x", result.x, err_x, "-");
	plt::named_plot("y", result.x, err_y, "-");
	plt::title("PID Controller error");
	plt::xlabel("t (s)");
	plt::ylabel("$Error\\:(\\circ)$");
	plt::legend();
	plt::xlim(0.0, result.x.back());

	plt::figure();
	plt::named_plot("x", result.x, theta_x, "-");
	plt::named_plot("y", result.x, theta_y, "-");
	plt::title("TVC");
	plt::xlabel("t (s)");
	plt::ylabel("$deflection\\:(\\circ)$");
	plt::legend();
	plt::xlim(0.0, result.x.back());

	plt::show();

	printf("Kp = [%.17e 0; 0 %.17e]\n", best[0], best[1]);
	printf("Ki = [%.17e 0; 0 %.17e]\n", best[2], best[3]);
	printf("Kd = [%.17e 0; 0 %.17e]\n", best[4], best[5]);

	std::cout << std::endl << "At 15 deg: ";
	objective(best, phi, theta, radians(15.0), tvc_bound, true);

	std::cout << std::endl << "At 235 deg: ";
	objective(best, phi, theta, radians(235.0), tvc_bound, true);

	std::cout << std::endl << "At 300 deg: ";
	objective(best, phi, theta, radians(300.0), tvc_bound, true);

	return 0;
}