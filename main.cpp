#include <iostream>

#include "math/gopt.hpp"

#include <chrono>
#include <iomanip>

#include "unit_tests.hpp"
#include <iostream>

using namespace gopt;

struct State
{
	union {
		struct {
			Vec3 dx;
			Vec3 omega;
			Quat q;
		};

		Vector<10> value;
	};

	State() {}
};

struct Observations
{
	union {
		struct {
			Vec3 accel;
			Vec3 gyro;
		};

		Vector<6> value;
	};

	Observations() {}
};

Vector<10> f(const Vector<10>& x, double t)
{
	State res;
	State& s = *(State*)(&x[0]);
	const gopt::Quat q = gopt::normalize(s.q);

	//////////////////////
	Vec2 tvc(0);
	struct {
		double m0 = 1.0;
		double d_com_ct = 0.45;
		double Ixx = 0.5;
		double Iyy = 0.5;
		double Izz = 0.2;
	} params;
	//////////////////////

	tvc[1] = radians(5.0) * std::sin(2 * gopt::pi<> * t);

	const double thrust_val = (t <= 5 ? 11 : 0);

	// Build thrust vector in local coordinates.
	const double sx = std::sin(tvc[0]);
	const double cx = std::cos(tvc[0]);
	const double sy = std::sin(tvc[1]);
	const double cy = std::cos(tvc[1]);
	const gopt::Vec3 thrust = gopt::Vec3(-sy, -sx * cy, cx * cy) * thrust_val;

	// Gravitational force vector
	const gopt::Vec3 weight = gopt::Vec3(0, 0, -params.m0 * 9.81);

	// Setup total forces and moments.
	const gopt::Vec3 forces = gopt::rotate(q, thrust) + weight;
	const gopt::Vec3 moments = gopt::cross(gopt::Vec3(0, 0, -params.d_com_ct), thrust);

	// Inertia matrix and its inverse.
	const gopt::Matrix<3, 3> I(params.Ixx, 0, 0, 0, params.Iyy, 0, 0, 0, params.Izz);
	const gopt::Matrix<3, 3> I_inv(1 / params.Ixx, 0, 0, 0, 1 / params.Iyy, 0, 0, 0, 1 / params.Izz);

	res.dx = forces / params.m0;
	res.omega = I_inv * (moments - gopt::cross(s.omega, I * s.omega));

	gopt::Quat omega;
	omega.w = 0;
	omega.v = s.omega;
	res.q = (omega * q) / 2;

	return res.value;
}

Vector<6> h(const Vector<10>& x, double t)
{
	Observations res;
	State& p = *(State*)(&x[0]);

	//////////////////////
	Vec2 tvc(0);
	struct {
		double m0 = 1.0;
	} params;
	//////////////////////

	tvc[1] = radians(5.0) * std::sin(2 * gopt::pi<> * t);

	const double thrust_val = (t <= 5 ? 11 : 0);

	// Build thrust vector in local coordinates.
	const double sx = std::sin(tvc[0]);
	const double cx = std::cos(tvc[0]);
	const double sy = std::sin(tvc[1]);
	const double cy = std::cos(tvc[1]);
	const gopt::Vec3 thrust = gopt::Vec3(-sy, -sx * cy, cx * cy) * thrust_val;

	// Gravitational force vector
	const gopt::Vec3 weight = gopt::Vec3(0, 0, -params.m0 * 9.81);

	// Setup total forces and moments.
	const gopt::Quat q = gopt::normalize(p.q);
	const gopt::Vec3 forces = gopt::rotate(q, thrust) + weight;
	res.accel = forces / params.m0;

	res.gyro = p.omega;

	return res.value;
}

int main()
{
#if _DEBUG
	unit_tests();
#endif
	State state;
	state.value = 0;
	state.q = Quaternion(1.0);

	Observations obs;
	obs.accel = { 0 };
	obs.gyro = { 0 };

	Matrix<10, 10> Pxx = gopt::eye<double, 10>() * 0.001;
	Matrix<10, 10> Q = gopt::eye<double, 10>() * 0.001;
	Matrix<6, 6> R = gopt::eye<double, 6>() * 0.001;

	auto start = std::chrono::system_clock::now();
	State cons;
	cons.value = state.value;
	const int N = 50000;
	const double h0 = 0.001;

	for (int i = 0; i < N; i++)
	{
		State curr_state;
		curr_state.value = state.value;

		curr_state.value = gopt::DP45(f, curr_state.value, i*h0, (i+1)*h0, 1e-6, 1e-6);
		obs.accel = (curr_state.dx - state.dx) / h0;
		obs.gyro = curr_state.omega;

		ekf(state.value, obs.value, Pxx, Q, R, f, h, i*h0, (i+1)*h0);
		state.q = gopt::normalize(state.q);
	}

	auto end = std::chrono::system_clock::now();
	const double time = static_cast<std::chrono::duration<double, std::nano>>(end - start).count() / 1e9;
	std::cout << "Time: " << time << " s, each: " << time / N << "s." << std::endl;
	std::cout << state.value << std::endl;
	std::cout << "Real: " << gopt::DP45(f, cons.value, 0.0, N * h0, 1e-9, 1e-9) << std::endl << std::endl;


	//const double time = static_cast<std::chrono::duration<double, std::nano>>(end - start).count() / 1e9;

	system("pause");
}