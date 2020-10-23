#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <string>
#include "constants.hpp"

namespace gopt
{
	template <typename T, unsigned int S>
	class Vector_t
	{
	private:
		std::array<T, S> data;

	private:
		Vector_t& add(const Vector_t& v)
		{
			for (int i = 0; i < S; i++)
				data[i] += v.data[i];
			return *this;
		}

		Vector_t& sub(const Vector_t& v)
		{
			for (int i = 0; i < S; i++)
				data[i] -= v.data[i];
			return *this;
		}

		T mul(const Vector_t& v) const
		{
			T res = 0;
			for (int i = 0; i < S; i++)
				res += data[i] * v.data[i];
			return res;
		}

		Vector_t& mul(const T& s)
		{
			for (int i = 0; i < S; i++)
				data[i] *= s;
			return *this;
		}

		Vector_t& div(const T& s)
		{
			for (int i = 0; i < S; i++)
				data[i] /= s;
			return *this;
		}

	public:
		Vector_t() {}

		Vector_t(const T& t)
		{
			data.fill(t);
		}

		template <typename... Ts, typename = typename std::enable_if_t<(sizeof...(Ts) == S)>>
		Vector_t(const Ts&... ts)
			: data{static_cast<T>(ts)...} {}

		template <unsigned int... Seq, typename Tuple>
		Vector_t(std::integer_sequence<unsigned int, Seq...>, Tuple&& tuple)
			: data{ static_cast<T>(std::get<Seq>(tuple))... } {}

		Vector_t(const Vector_t& v)
			: data(v.data) {}

		template <typename V, typename = std::enable_if_t<!std::is_same_v<T, V>>>
		Vector_t(const Vector_t<V, S>& v)
		{
			for (int i = 0; i < S; i++)
				data[i] = static_cast<T>(v[i]);
		}

		friend Vector_t operator-(const Vector_t& v)
		{
			Vector_t<T, S> res;
			for (int i = 0; i < S; i++)
				res.data[i] = -v.data[i];
			return res;
		}

		friend Vector_t operator+(Vector_t rhs, const Vector_t& lhs)
		{
			return rhs.add(lhs);
		}

		friend Vector_t operator-(Vector_t rhs, const Vector_t& lhs)
		{
			return rhs.sub(lhs);
		}

		friend T operator*(const Vector_t& lhs, const Vector_t& rhs)
		{
			return lhs.mul(rhs);
		}

		friend Vector_t operator*(const T& s, Vector_t v)
		{
			return v.mul(s);
		}

		friend Vector_t operator*(Vector_t v, const T& s)
		{
			return v.mul(s);
		}

		friend Vector_t operator/(Vector_t v, const T& s)
		{
			return v.div(s);
		}

		Vector_t& operator+=(const Vector_t& v)
		{
			return add(v);
		}

		Vector_t& operator-=(const Vector_t& v)
		{
			return sub(v);
		}

		Vector_t& operator*=(const T& s)
		{
			return add(s);
		}

		Vector_t& operator/=(const T& s)
		{
			return div(s);
		}

		bool operator==(const Vector_t& v) const
		{
			for (int i = 0; i < S; i++)
				if (data[i] != v.data[i])
					return false;
			return true;
		}

		bool operator!=(const Vector_t& v) const
		{
			for (int i = 0; i < S; i++)
				if (data[i] != v.data[i])
					return true;
			return false;
		}

		bool almost_equal(const Vector_t& v) const
		{
			for (int i = 0; i < S; i++)
				if (std::abs(data[i] - v.data[i]) > weak_epsilon<T>)
					return false;
			return true;
		}

		T magnitude() const
		{
			T res = 0;
			for (int i = 0; i < S; i++)
				res += data[i] * data[i];
			return res;
		}

		T length() const
		{
			return std::sqrt(magnitude());
		}

		T* begin()
		{
			return &(data[0]);
		}

		const T* begin() const
		{
			return &(data[0]);
		}

		T* end()
		{
			return &(data[0]) + S;
		}

		const T* end() const
		{
			return &(data[0]) + S;
		}

		T& operator[](const unsigned int index)
		{
			return data[index];
		}

		const T& operator[](const unsigned int index) const
		{
			return data[index];
		}

		std::string toString() const
		{
			std::string res = "Vector<";
			res += std::string(typeid(T).name()) + "," + std::to_string(S) + ">{" + std::to_string(data[0]);
			for (int i = 1; i < S; i++)
				res += ", " + std::to_string(data[i]);
			res += "}";
			return res;
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector_t<T, S>& v)
		{
			os << v.toString();
			return os;
		}

		constexpr unsigned int size() const { return S; }
	};

#ifdef SINGLE_PRECISION
	template <unsigned int S>
	using Vector = Vector_t<float, S>;
#else
	template <unsigned int S>
	using Vector = Vector_t<double, S>;
#endif

	using Vec2 = Vector_t<double, 2>;
	using Vec3 = Vector_t<double, 3>;
	using Vec4 = Vector_t<double, 4>;

	using Vec2f = Vector_t<float, 2>;
	using Vec3f = Vector_t<float, 3>;
	using Vec4f = Vector_t<float, 4>;
}