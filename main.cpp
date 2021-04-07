#include <iostream>
#include "math/gopt.hpp"
using namespace gopt;

double gravity(const double h)
{
	const double Rt = 6371000;
	const double g0 = 9.81;
	return g0 * std::pow(Rt / (Rt + h), 2);
}

double density_ISA(const double h)
{
	const double g = gravity(h);
	const double R = 286.9;

	if (h <= 11000)
	{
		const double rho0 = 1.225;
		const double T0 = 288.15;
		const double alpha = -0.0065;

		const double T = T0 + alpha * h;
		return rho0 * std::pow(T / T0, -1 - g / (alpha * R));
	}
	else if (h <= 20000)
	{
		const double h0 = 11000;
		const double rho0 = 0.365324;
		const double T0 = 216.65;

		return rho0 * std::exp(-g / (R * T0) * (h - h0));
	}
	else if (h <= 32000)
	{
		const double h0 = 20000;
		const double rho0 = 0.0890522;
		const double T0 = 216.65;
		const double alpha = 0.001;

		const double T = T0 + alpha * (h - h0);
		return rho0 * std::pow(T / T0, -1 - g / (alpha * R));
	}
	else if (h <= 47000)
	{
		const double h0 = 32000;
		const double rho0 = 0.013604;
		const double T0 = 228.65;
		const double alpha = 0.0028;

		const double T = T0 + alpha * (h - h0);
		return rho0 * std::pow(T / T0, -1 - g / (alpha * R));
	}
	else if (h <= 51000)
	{
		const double h0 = 47000;
		const double rho0 = 0.00151052;
		const double T0 = 270.65;

		return rho0 * std::exp(-g / (R * T0) * (h - h0));
	}
	else if (h <= 71000)
	{
		const double h0 = 51000;
		const double rho0 = 0.000918605;
		const double T0 = 270.65;
		const double alpha = -0.0028;

		const double T = T0 + alpha * (h - h0);
		return rho0 * std::pow(T / T0, -1 - g / (alpha * R));
	}
	else if (h <= 84852)
	{
		const double h0 = 71000;
		const double rho0 = 7.26618e-05;
		const double T0 = 214.65;
		const double alpha = -0.002;

		const double T = T0 + alpha * (h - h0);
		return rho0 * std::pow(T / T0, -1 - g / (alpha * R));
	}

	return 0.0;
}

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
	Vec2 tvc0;
	double tvc_dot;
	Vec2 tvc_dot_sign;
	double t0;
	Vec2 period;
	double tb;
	double d_com_cp;
	double Cd0;
	double Cda;
	double Cla;
	double rho;
	double Sref;
	double T;
};

Vec2 calculate_tvc(const double& t, const Params& p)
{
	const double dt = t - p.t0;
	Vec2 res;

	res[0] = (dt < p.period[0]) ? (p.tvc0[0] + p.tvc_dot_sign[0] * dt) : p.tvc[0];
	res[1] = (dt < p.period[1]) ? (p.tvc0[1] + p.tvc_dot_sign[1] * dt) : p.tvc[1];
	
	return res;
}

Vector_t<double, 13> dynamics(const Vector_t<double, 13>& x, const double& t, const Params& p)
{
	State& s = *(State*)(&x[0]);

	const double T_val = (t < p.tb) ? p.T : 0;

	const Vec2 tvc = calculate_tvc(t, p);

	const double sx = std::sin(tvc[0]);
	const double cx = std::cos(tvc[0]);
	const double sy = std::sin(tvc[1]);
	const double cy = std::cos(tvc[1]);
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

		const double rho = density_ISA(s.x[2]);
		const Vec3 L = 0.5 * rho * v_mag * v_mag * p.Cla * alpha * p.Sref * nL;
		const Vec3 D = 0.5 * rho * v_mag * v_mag * (p.Cd0 + p.Cda * alpha) * p.Sref * nD;
		F_ae = L + D;
	}
	else
		F_ae = Vec3(0);

	const Vec3 forces = rotate(conjugate(s.q), T + F_ae) + Vec3(0, 0, -m * gravity(s.x[2]));
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

	Quat qx = Quat(radians(10.0), Vec3(1, 0, 0));
	Quat qy = Quat(radians(10.0), Vec3(0, 1, 0));
	Quat qz = Quat(radians(10.0), Vec3(0, 0, 1));
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
	params.tvc0 = 0; // rad
	params.tvc_dot = radians(500.0); // rad/s
	params.tvc_dot_sign = 0;
	params.t0 = 0; // s
	params.period = 0; // s
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
		params.t0 = prev_t;
		const Vec2 period = (params.tvc - params.tvc0) / params.tvc_dot;
		params.period = abs(period);
		params.tvc_dot_sign = params.tvc_dot;
		if (period[0] < 0)
			params.tvc_dot_sign[0] = -params.tvc_dot;
		if (period[1] < 0)
			params.tvc_dot_sign[1] = -params.tvc_dot;

		t = prev_t + dt;
		const auto fcn = [&](const Vector_t<double, 13> x, double t) { return dynamics(x, t, params); };
		const auto res = DP45(fcn, state.value, prev_t, t, 1e-9, 1e-9);

		// Finished
		params.tvc0 = calculate_tvc(t, params);

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

			v[15] = params.tvc0[0];
			v[16] = params.tvc0[1];

			// Dynamic pressure
			v[17] = 0.5 * density_ISA(state.x[2]) * state.v.magnitude();

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