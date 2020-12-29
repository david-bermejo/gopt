#pragma once

#include <cmath>
#include <string>
#include <iostream>

#include "algorithms.hpp"
#include "vector.hpp"

namespace gopt
{
	template <typename T>
	class Quaternion_t
	{
	public:
		T w;
		union {
			struct { T x, y, z; };
			gopt::Vector_t<T, 3> v;
		};

	private:
		Quaternion_t& add(const Quaternion_t& q)
		{
			w += q.w;
			v += q.v;

			return *this;
		}

		Quaternion_t& sub(const Quaternion_t& q)
		{
			w -= q.w;
			v -= q.v;

			return *this;
		}

		Quaternion_t& mul(const Quaternion_t& q)
		{
			Quaternion_t res;

			res.w = w * q.w - dot(v, q.v);
			res.v = w * q.v + q.w * v + cross(v, q.v);

			w = res.w;
			v = res.v;

			return *this;
		}

		Quaternion_t& mul(const T& s)
		{
			w *= s;
			v *= s;

			return *this;
		}

		Quaternion_t& div(const Quaternion_t& q)
		{
			return mul(inverse(q));
		}

		Quaternion_t& div(const T& s)
		{
			w /= s;
			v /= s;

			return *this;
		}

	public:
		Quaternion_t() {}

		Quaternion_t(const T& s) : w(s), v(static_cast<T>(0)) {}

		template <typename V, typename... Ts, typename = typename std::enable_if_t<std::conjunction_v<std::is_same<V, Ts>...> && (sizeof...(Ts) == 3)>>
		Quaternion_t(const V& first, const Ts&... ts)
			: w(static_cast<T>(first)), v{ ts... } {}

		// Positive angle produces counter clockwise rotation around axis (right-hand rule).
		Quaternion_t(const T& angle, const Vector_t<T, 3>& axis)
			: w(std::cos(angle/2)), v(axis * std::sin(-angle/2)) {}

		Quaternion_t(const Quaternion_t& q)
			: w(q.w), v(q.v) {}

		template <typename V, typename = std::enable_if_t<!std::is_same_v<T, V>>>
		Quaternion_t(const Quaternion_t<V>& q)
			: w(q.w), v(q.v) {}

		friend Quaternion_t operator-(const Quaternion_t& q)
		{
			return Quaternion_t(-q.w, -q.x, -q.y, -q.z);
		}

		friend Quaternion_t operator+(Quaternion_t lhs, const Quaternion_t& rhs)
		{
			return lhs.add(rhs);
		}

		friend Quaternion_t operator-(Quaternion_t lhs, const Quaternion_t& rhs)
		{
			return lhs.sub(rhs);
		}

		friend Quaternion_t operator*(Quaternion_t lhs, const Quaternion_t& rhs)
		{
			return lhs.mul(rhs);
		}

		friend Quaternion_t operator*(Quaternion_t lhs, const T& rhs)
		{
			return lhs.mul(rhs);
		}

		friend Quaternion_t operator/(Quaternion_t lhs, const Quaternion_t& rhs)
		{
			return lhs.div(rhs);
		}

		friend Quaternion_t operator/(Quaternion_t lhs, const T& rhs)
		{
			return lhs.div(rhs);
		}

		Quaternion_t& operator+=(const Quaternion_t& q)
		{
			return add(q);
		}

		Quaternion_t& operator-=(const Quaternion_t& q)
		{
			return sub(q);
		}

		Quaternion_t& operator*=(const Quaternion_t& q)
		{
			return mul(q);
		}

		Quaternion_t& operator*=(const T& s)
		{
			return mul(s);
		}

		Quaternion_t& operator/=(const Quaternion_t& q)
		{
			return div(q);
		}

		Quaternion_t& operator/=(const T& s)
		{
			return div(s);
		}

		bool operator==(const Quaternion_t& q) const
		{
			if (w != q.w)
				return false;

			for (int i = 0; i < 3; i++)
				if (v[i] != q.v[i])
					return false;

			return true;
		}

		bool operator!=(const Quaternion_t& q) const
		{
			if (w != q.w)
				return true;

			for (int i = 0; i < 3; i++)
				if (v[i] != q.v[i])
					return true;

			return false;
		}

		bool almost_equal(const Quaternion_t& q) const
		{
			if (std::abs(w - q.w) > weak_epsilon<T>)
				return false;

			return v.almost_equal(q.v);
		}

		inline T magnitude() const
		{
			return w*w + v.magnitude();
		}

		inline T length() const
		{
			return std::sqrt(magnitude());
		}

		T& operator[](const unsigned int index)
		{
			return *(T*)(&w + index);
		}

		const T& operator[](const unsigned int index) const
		{
			return *(T*)(&w + index);
		}

		std::string toString() const
		{
			std::string res = "Quaternion<";
			res += std::string(typeid(T).name()) + ">{" + std::to_string(w);
			for (int i = 0; i < 3; i++)
				res += ", " + std::to_string(v[i]);
			res += "}";
			return res;
		}

		friend std::ostream& operator<<(std::ostream& os, const Quaternion_t& q)
		{
			os << q.toString();
			return os;
		}

		constexpr unsigned int size() const { return 4; }
	};

	#ifdef SINGLE_PRECISION
	using Quaternion = Quaternion_t<float>;
	#else
	using Quaternion = Quaternion_t<double>;
	#endif

	using Quat = Quaternion;
};