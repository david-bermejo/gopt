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
		friend Matrix_t<T, R, Cols> operator*(const Matrix_t<T, R, C>& lhs, const Matrix_t<T, C, Cols>& rhs)
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
			return &data[0][0] + R * C;
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
	using Mat = Matrix_t<float, R, C>;
#else
	template <unsigned int R, unsigned int C>
	using Mat = Matrix_t<double, R, C>;
#endif

	using Mat2 = Matrix_t<double, 2, 2>;
	using Mat3 = Matrix_t<double, 3, 3>;
	using Mat4 = Matrix_t<double, 4, 4>;

	using Mat2f = Matrix_t<float, 2, 2>;
	using Mat3f = Matrix_t<float, 3, 3>;
	using Mat4f = Matrix_t<float, 4, 4>;

	template <typename T>
	class Matrix
	{
	private:
		unsigned int _rows = 0;
		unsigned int _cols = 0;
		T** data = nullptr;

	private:
		Matrix& add(const Matrix& m)
		{
			assert(_rows == m._rows && _cols == m._cols);

			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					data[i][j] += m.data[i][j];

			return *this;
		}

		Matrix& sub(const Matrix& m)
		{
			assert(_rows == m._rows && _cols == m._cols);

			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					data[i][j] -= m.data[i][j];
			
			return *this;
		}

		Matrix& mul(const Matrix& m)
		{
			assert(_cols == m._rows);

			Matrix<T> res(_rows, m._cols);

			for (int i = 0; i < _rows; i++)
			{
				for (int j = 0; j < m._cols; j++)
				{
					T tmp = 0;
					for (int k = 0; k < _cols; k++)
						tmp += data[i][k] * m.data[k][j];
					
					res.data[i][j] = tmp;
				}
			}

			return move(res);
		}

		Matrix& mul(const T& s)
		{
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					data[i][j] *= s;
			
			return *this;
		}

		Matrix& div(const T& s)
		{
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					data[i][j] /= s;
			
			return *this;
		}

	public:
		Matrix() {}

		Matrix(unsigned int rows, unsigned int cols)
			: _rows(rows), _cols(cols), data(new T*[rows])
		{
			data[0] = new T[rows * cols];

			T* const beg = data[0];
			for (int i = 1; i < rows; i++)
				data[i] = beg + i * cols;
		}

		Matrix(const Matrix& m)
			: _rows(m._rows), _cols(m._cols), data(new T*[m._rows])
		{
			data[0] = new T[_rows * _cols];
			std::copy(m.data[0], &m.data[0][_rows * _cols], data[0]);

			T* const beg = data[0];
			for (int i = 1; i < _rows; i++)
				data[i] = beg + i * _cols;
		}

		Matrix(Matrix&& m)
			: _rows(m._rows), _cols(m._cols), data(m.data)
		{
			m._rows = 0;
			m._cols = 0;
			m.data = nullptr;
		}

		~Matrix()
		{
			if (data)
			{
				delete[] data[0];
				delete[] data;
			}
		}

		Matrix& operator=(const Matrix& m)
		{
			_rows = m._rows;
			_cols = m._cols;

			if (data) {
				delete[] data[0];
				delete[] data;
			}
			
			data = new T*[_rows];
			data[0] = new T[_rows * _cols];
			std::copy(m.data[0], &m.data[0][_rows * _cols], data[0]);

			T* const beg = data[0];
			for (int i = 1; i < _rows; i++)
				data[i] = beg + i * _cols;
			
			return *this;
		}

		Matrix& fill(const T& s)
		{
			std::fill_n(data[0], _rows * _cols, s);
			return *this;
		}

		template <typename... Ts, typename = std::enable_if_t<(sizeof...(Ts) > 1)>>
		Matrix& fill(Ts...  ts)
		{
			assert(_rows * _cols == sizeof...(ts));

			unsigned int index = 0;
			((data[0][index++] = static_cast<T>(ts)), ...);

			return *this;
		}

		Matrix& fill(const std::vector<T>& src)
		{
			std::copy(src.data(), src.data() + _rows * _cols, data[0]);
			return *this;
		}

		template <typename F, typename... Args, typename = std::enable_if_t<std::is_invocable_v<F, Args...>>>
		Matrix& fill(F&& f, Args... args)
		{
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					data[i][j] = f(args...);
			return *this;
		}

		Matrix& move(Matrix& m)
		{
			_rows = m._rows;
			_cols = m._cols;

			if (data)
			{
				delete[] data[0];
				delete[] data;
			}

			data = m.data;
			m.data = nullptr;

			m._rows = 0;
			m._cols = 0;

			return *this;
		}

		T& operator()(unsigned int i, unsigned int j)
		{
			assert(i < _rows && j < _cols);
			return data[i][j];
		}

		const T& operator()(unsigned int i, unsigned int j) const
		{
			assert(i < _rows && j < _cols);
			return data[i][j];
		}

		Vector<T> operator()(unsigned int i) const
		{
			assert(i < _rows);

			Vector<T> res;
			res.len = _cols;
			res._destroyable = false;
			res.data = data[i];
			
			return res;
		}

		friend Matrix operator-(const Matrix& m)
		{
			Matrix<T> res(m._rows, m._cols);

			for (int i = 0; i < m._rows; i++)
				for (int j = 0; j < m._cols; j++)
					res.data[i][j] = -m.data[i][j];
			return res;
		}

		friend Matrix operator+(Matrix lhs, const Matrix& rhs)
		{
			return lhs.add(rhs);
		}

		friend Matrix operator-(Matrix lhs, const Matrix& rhs)
		{
			return lhs.sub(rhs);
		}

		friend Matrix operator*(const Matrix& lhs, const Matrix& rhs)
		{
			assert(lhs._cols == rhs._rows);

			Matrix<T> res(lhs._rows, rhs._cols);

			for (int i = 0; i < lhs._rows; i++)
			{
				for (int j = 0; j < rhs._cols; j++)
				{
					T tmp = 0;
					for (int k = 0; k < lhs._cols; k++)
						tmp += lhs.data[i][k] * rhs.data[k][j];
					
					res.data[i][j] = tmp;
				}
			}

			return res;
		}

		friend Matrix operator*(const T& s, Matrix m)
		{
			return m.mul(s);
		}

		friend Matrix operator*(Matrix m, const T& s)
		{
			return m.mul(s);
		}

		friend Matrix operator/(const T& s, Matrix m)
		{
			return m.div(s);
		}

		friend Matrix operator/(Matrix m, const T& s)
		{
			return m.div(s);
		}

		Matrix& operator+=(const Matrix& m)
		{
			return add(m);
		}

		Matrix& operator-=(const Matrix& m)
		{
			return sub(m);
		}

		Matrix& operator*=(const Matrix& m)
		{
			return mul(m);
		}

		Matrix& operator*=(const T& s)
		{
			return mul(s);
		}

		Matrix& operator/=(const T& s)
		{
			return div(s);
		}

		bool operator==(const Matrix& m) const
		{
			if (_rows != m._rows || _cols != m._cols)
				return false;
			
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					if (data[i][j] != m.data[i][j])
						return false;
			return true;
		}

		bool operator!=(const Matrix& m) const
		{
			if (_rows != m._rows || _cols != m._cols)
				return true;
			
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					if (data[i][j] != m.data[i][j])
						return true;
			return false;
		}

		bool almost_equal(const Matrix& m) const
		{
			if (_rows != m._rows || _cols != m._cols)
				return false;

			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					if (std::abs(data[i][j] - m.data[i][j]) > weak_epsilon<T>)
						return false;
			return true;
		}

		T* begin()
		{
			return data[0];
		}

		const T* begin() const
		{
			return data[0];
		}

		T* end()
		{
			return data[0] + _rows * _cols;
		}

		const T* end() const
		{
			return data[0] * _rows * _cols;
		}

		T* operator[](unsigned int index)
		{
			assert(index < _rows);
			return data[index];
		}

		const T* operator[](unsigned int index) const
		{
			assert(index < _rows);
			return data[index];
		}

		std::string toString() const
		{
			std::string res = "Matrix<";
			res += std::string(typeid(T).name()) + ", " + std::to_string(_rows) + ", " + std::to_string(_cols) + ">{";
			for (int i = 0; i < _rows; i++)
			{
				res += "\n\t[" + std::to_string(data[i][0]);

				for (int j = 1; j < _cols; j++)
					res += ", " + std::to_string(data[i][j]);

				res += "]";
			}
			res += "\n}";

			return res;
		}

		friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& v)
		{
			os << v.toString();
			return os;
		}

		const unsigned int rows() const { return _rows; }
		const unsigned int cols() const { return _cols; }
		const unsigned int size() const { return _rows * _cols; }
	};
}