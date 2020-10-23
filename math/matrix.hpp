#pragma once

#include "vector.hpp"
#include "../meta/range.hpp"

#include <algorithm>
#include <string>

namespace gopt
{
	template <typename T, unsigned int R, unsigned int C>
	class Matrix_t
	{
	private:
		Vector_t<Vector_t<T, C>, R> data;

	private:
		template <unsigned int... Seq, typename Tuple>
		Matrix_t(std::integer_sequence<unsigned int, Seq...>, Tuple&& tuple)
			: data{Vector_t<T, C>(meta::make_integer_range<Seq * C, (Seq + 1) * C>{}, tuple)...} {}
		
		Matrix_t& add(const Matrix_t& m)
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					data[i][j] += m.data[i][j];
			return *this;
		}

		Matrix_t& sub(const Matrix_t& m)
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					data[i][j] -= m.data[i][j];
			return *this;
		}

		template <typename = typename std::enable_if_t<R == C>>
		Matrix_t& mul(const Matrix_t& m)
		{
			Matrix_t<T, R, C> res;

			for (int i = 0; i < R; i++)
			{
				for (int j = 0; j < C; j++)
				{
					T tmp = 0;
					for (int k = 0; k < R; k++)
						tmp += data[i][k] * m.data[k][j];
					res[i][j] = tmp;
				}
			}

			data = res.data;
			return *this;
		}

		Matrix_t& mul(const T& s)
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					data[i][j] *= s;
			return *this;
		}

		Matrix_t& div(const T& s)
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					data[i][j] /= s;
			return *this;
		}

	public:
		Matrix_t() {}

		Matrix_t(const T& t)
		{
			std::fill_n(&data[0][0], R * C, t);
		}

		template <typename... Ts, typename = typename std::enable_if_t<(sizeof...(Ts) == R*C)>>
		Matrix_t(Ts... ts)
			: Matrix_t(std::make_integer_sequence<unsigned int, R>{}, std::make_tuple(ts...)) {}

		Matrix_t(const Matrix_t& m)
			: data(m.data) {}

		template <typename V, typename = std::enable_if_t<!std::is_same_v<T, V>>>
		Matrix_t(const Matrix_t<V, R, C>& m)
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					data[i][j] = static_cast<T>(m[i][j]);
		}

		friend Matrix_t operator+(Matrix_t lhs, const Matrix_t& rhs)
		{
			return lhs.add(rhs);
		}

		friend Matrix_t operator-(Matrix_t lhs, const Matrix_t& rhs)
		{
			return lhs.sub(rhs);
		}

		template <unsigned int Cols>
		friend Matrix_t operator*(const Matrix_t<T, R, C>& lhs, const Matrix_t<T, C, Cols>& rhs)
		{
			Matrix_t<T, R, Cols> res;

			for (int i = 0; i < R; i++)
			{
				for (int j = 0; j < Cols; j++)
				{
					T tmp = 0;
					for (int k = 0; k < C; k++)
						tmp += lhs[i][k] * rhs[k][j];
					res[i][j] = tmp;
				}
			}

			return res;
		}

		friend Matrix_t operator*(const T& s, Matrix_t m)
		{
			return m.mul(s);
		}

		friend Matrix_t operator*(Matrix_t m, const T& s)
		{
			return m.mul(s);
		}

		friend Matrix_t operator/(Matrix_t m, const T& s)
		{
			return m.div(s);
		}

		Matrix_t& operator+=(const Matrix_t& m)
		{
			return add(m);
		}

		Matrix_t& operator-=(const Matrix_t& m)
		{
			return sub(m);
		}

		Matrix_t& operator*=(const Matrix_t& m)
		{
			return mul(m);
		}

		Matrix_t& operator*=(const T& s)
		{
			return mul(s);
		}

		Matrix_t& operator/=(const T& s)
		{
			return div(s);
		}

		bool operator==(const Matrix_t& m) const
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					if (data[i][j] != m.data[i][j])
						return false;
			return true;
		}

		bool operator!=(const Matrix_t& m) const
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					if (data[i][j] != m.data[i][j])
						return true;
			return false;
		}

		bool almost_equal(const Matrix_t& m) const
		{
			for (int i = 0; i < R; i++)
				for (int j = 0; j < C; j++)
					if (std::abs(data[i][j] - m.data[i][j]) > weak_epsilon<T>)
						return false;
			return true;
		}

		T* begin()
		{
			return &data[0][0];
		}

		const T* begin() const
		{
			return &data[0][0];
		}

		T* end()
		{
			return &data[0][0] + R*C;
		}

		const T* end() const
		{
			return &data[0][0] + R * C;
		}

		Vector_t<T, C>& operator[](unsigned int index)
		{
			return data[index];
		}

		const Vector_t<T, C>& operator[](unsigned int index) const
		{
			return data[index];
		}

		std::string toString() const
		{
			std::string res = "Matrix<";
			res += std::string(typeid(T).name()) + ", " + std::to_string(R) + ", " + std::to_string(C) + ">{";
			for (int i = 0; i < R; i++)
			{
				res += "\n\t[" + std::to_string(data[i][0]);

				for (int j = 1; j < C; j++)
					res += ", " + std::to_string(data[i][j]);

				res += "]";
			}
			res += "\n}";

			return res;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix_t<T, R, C>& v)
		{
			os << v.toString();
			return os;
		}

		constexpr unsigned int rows() const { return R; }
		constexpr unsigned int columns() const { return C; }
		constexpr unsigned int size() const { return R * C; }
	};

#ifdef SINGLE_PRECISION
	template <unsigned int R, unsigned int C>
	using Matrix = Matrix_t<float, R, C>;
#else
	template <unsigned int R, unsigned int C>
	using Matrix = Matrix_t<double, R, C>;
#endif

	using Mat3 = Matrix_t<double, 3, 3>;
	using Mat4 = Matrix_t<double, 4, 4>;

	using Mat3f = Matrix_t<float, 3, 3>;
	using Mat4f = Matrix_t<float, 4, 4>;
}