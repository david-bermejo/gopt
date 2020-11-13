#pragma once

namespace gopt
{
	template <typename T>
	class PID
	{
	private:
		T Kp, Ki, Kd;
		T set_point;

		T prev_deriv;
		T int_accum;
		T prev_err;

	public:
		PID() = default;

		PID(const T& Kp, const T& Ki, const T& Kd)
			: Kp(Kp), Ki(Ki), Kd(Kd), set_point(0), prev_deriv(0), int_accum(0), prev_err(0) {}

		void set(const T& sp)
		{
			set_point = sp;
			prev_err = sp;
		}

		T update(const T& input, const T dt)
		{
			T output;

			// Error = sp - input
			const T err = set_point - input;

			// Proportional term
			output = Kp * err;

			// Integral term (simple trapezoidal formula)
			int_accum += dt/2 * (err + prev_err);
			output += Ki * int_accum;

			// Derivative term
			output += Kd * (err - prev_err) / dt;
			prev_err = err;

			return output;
		}
	};
}