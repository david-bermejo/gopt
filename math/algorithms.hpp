#pragma once

#include "vector.hpp"
#include "matrix.hpp"
#include "constants.hpp"

namespace gopt
{
	template <typename T, unsigned int S>
	T max(const Vector_t<T, S>& v)
	{
		T res = v[0];

		for (int i = 1; i < S; i++)
			if (res < v[i])
				res = v[i];

		return res;
	}

	template <typename T, unsigned int S>
	T min(const Vector_t<T, S>& v)
	{
		T res = v[0];

		for (int i = 1; i < S; i++)
			if (res > v[i])
				res = v[i];

		return res;
	}

	template <typename T, unsigned int S>
	T dot(const Vector_t<T, S>& lhs, const Vector_t<T, S>& rhs)
	{
		return lhs * rhs;
	}

	template <typename T, unsigned int S>
	Matrix_t<T, S, S> outer(const Vector_t<T, S>& lhs, const Vector_t<T, S>& rhs)
	{
		Matrix_t<T, S, S> res;

		for (int i = 0; i < S; i++)
			for (int j = 0; j < S; j++)
				res[i][j] = lhs[i] * rhs[j];

		return res;
	}

	template <typename T>
	Vector_t<T, 3> cross(const Vector_t<T, 3>& lhs, const Vector_t<T, 3>& rhs)
	{
		return Vector_t<T, 3>
		{
			lhs[1] * rhs[2] - lhs[2] * rhs[1],
			lhs[2] * rhs[0] - lhs[0] * rhs[2],
			lhs[0] * rhs[1] - lhs[1] * rhs[0]
		};
	}

	template <typename T, unsigned int N>
	Vector_t<T, N> abs(Vector_t<T, N> v)
	{
		for (auto& i : v)
			i = std::abs(i);

		return v;
	}

	template <typename T, unsigned int N>
	unsigned int maxloc(const Vector_t<T, N>& v)
	{
		T max = v[0];
		int index = 0;

		for (int i = 1; i < N; i++)
		{
			if (v[i] > max)
			{
				index = i;
				max = v[i];
			}
		}

		return index;
	}

	template <typename T, unsigned int S>
	T angle(const Vector_t<T, S>& lhs, const Vector_t<T, S>& rhs)
	{
		return std::acos(dot(lhs, rhs) / (lhs.length() * rhs.length()));
	}

	template <typename T, unsigned int R, unsigned int C>
	Vector_t<T, R> operator*(const Matrix_t<T, R, C>& m, const Vector_t<T, C>& v)
	{
		Vector_t<T, R> res;

		for (int i = 0; i < R; i++)
		{
			T tmp = 0;
			for (int j = 0; j < C; j++)
				tmp += m[i][j] * v[j];
			res[i] = tmp;
		}

		return res;
	}

	template <typename T, unsigned int R, unsigned int C>
	Vector_t<T, C> operator*(const Vector_t<T, R>& v, const Matrix_t<T, R, C>& m)
	{
		Vector_t<T, C> res;

		for (int i = 0; i < C; i++)
		{
			T tmp = 0;
			for (int j = 0; j < R; j++)
				tmp += v[j] * m[j][i];
			res[i] = tmp;
		}

		return res;
	}

	template <typename T, unsigned int N>
	Matrix_t<T, N, N> eye()
	{
		Matrix_t<T, N, N> res(0);

		for (int i = 0; i < N; i++)
			res[i][i] = 1;

		return res;
	}

	template <typename T>
	Matrix_t<T, 4, 4> perspective(T ar, T fov, T z_near, T z_far)
	{
		const T tan_fov = std::tan(radians(fov / 2));
		const T sub = z_near - z_far;
		
		Matrix_t<T, 4, 4> res(0);
		res[0][0] = static_cast<T>(1) / (ar * tan_fov);
		res[1][1] = static_cast<T>(1) / tan_fov;
		res[2][2] = -(z_near + z_far) / sub;
		res[2][3] = -2 * z_near * z_far / sub;
		res[3][2] = -1;

		return res;
	}

	template <typename T>
	Matrix_t<T, 4, 4> translation(const T& x, const T& y, const T& z)
	{
		Matrix_t<T, 4, 4> res = eye<T, 4>();
		res[0][3] = x;
		res[1][3] = y;
		res[2][3] = z;

		return res;
	}

	template <typename T>
	Matrix_t<T, 4, 4> translation(const Vector_t<T, 3>& v)
	{
		Matrix_t<T, 4, 4> res = eye<T, 4>();
		res[0][3] = v[0];
		res[1][3] = v[1];
		res[2][3] = v[2];

		return res;
	}

	template <typename T>
	Matrix_t<T, 4, 4> scaling(const T& x, const T& y, const T& z)
	{
		Matrix_t<T, 4, 4> res = eye<T, 4>();
		res[0][0] = x;
		res[1][1] = y;
		res[2][2] = z;

		return res;
	}

	template <typename T>
	Matrix_t<T, 4, 4> scaling(const Vector_t<T, 3>& v)
	{
		Matrix_t<T, 4, 4> res = eye<T, 4>();
		res[0][0] = v[0];
		res[1][1] = v[1];
		res[2][2] = v[2];

		return res;
	}

	// Build ZYX rotation matrix using yaw, pitch and roll angles (Euler angles).
	template <typename T>
	Matrix_t<T, 4, 4> rotation(const T& yaw, const T& pitch, const T& roll)
	{
		const T sz = std::sin(yaw);
		const T cz = std::cos(yaw);
		const T sy = std::sin(pitch);
		const T cy = std::cos(pitch);
		const T sx = std::sin(roll);
		const T cx = std::cos(roll);

		return Matrix_t<T, 4, 4>
		{
			cz*cy, cz*sy*sx - sz*cx, cz*sy*cx + sz*sx, 0,
			sz*cy, sz*sy*sx + cz*cx, sz*sy*cx - cz*sx, 0,
			-sy, cy*sx, cy*cx, 0,
			0, 0, 0, 1
		};
	}

	// Return Euler angles from a ZYX rotation matrix in the form of Vec3{yaw, pitch, roll}.
	template <typename T>
	Vector_t<T, 3> euler_angles(const Matrix_t<T, 4, 4>& m)
	{
		const T r11 = m[0][0];
		const T r21 = m[1][1];

		return Vector_t<T, 3>
		{
			std::atan2(r21, r11),
			std::atan2(-m[2][0], std::sqrt(r11*r11 + r21*r21)),
			std::atan2(m[2][1], m[2][2])
		};
	}

	template <typename T, unsigned int R, unsigned int C>
	Matrix_t<T, C, R> transpose(const Matrix_t<T, R, C>& m)
	{
		Matrix_t<T, C, R> res;

		for (int i = 0; i < R; i++)
			for (int j = 0; j < C; j++)
				res[j][i] = m[i][j];

		return res;
	}

	template <typename T, unsigned int R, unsigned int C>
	void swap_row(Matrix_t<T, R, C>& M, unsigned int a, unsigned int b)
	{
		for (int i = 0; i < C; i++)
		{
			T tmp = M[a][i];
			M[a][i] = M[b][i];
			M[b][i] = tmp;
		}
	}

	template <typename T, unsigned int R, unsigned int C>
	void swap_column(Matrix_t<T, R, C>& M, unsigned int a, unsigned int b)
	{
		for (int i = 0; i < R; i++)
		{
			T tmp = M[i][a];
			M[i][a] = M[i][b];
			M[i][b] = tmp;
		}
	}

	template <typename T, unsigned int R, unsigned int C>
	Vector_t<T, C> row(const Matrix_t<T, R, C>& M, unsigned int index)
	{
		return M[index];
	}

	template <typename T, unsigned int R, unsigned int C>
	Vector_t<T, R> column(const Matrix_t<T, R, C>& M, unsigned int index)
	{
		Vector_t<T, R> res;

		for (int i = 0; i < R; i++)
			res[i] = M[i][index];

		return res;
	}

	template <typename T, unsigned int R, unsigned int C>
	void lu(const Matrix_t<T, R, C>& M, Matrix_t<T, R, C>& L, Matrix_t<T, R, C>& U)
	{
		for (int i = 0; i < R; i++)
		{
			L[i][i] = 1;

			for (int j = 0; j < i; j++)
				U[i][j] = 0;

			for (int j = i+1; j < C; j++)
				L[i][j] = 0;

			for (int k = i; k < C; k++)
			{
				T sum = 0;
				for (int j = 0; j < i; j++)
					sum += L[i][j] * U[j][k];
				U[i][k] = M[i][k] - sum;
			}

			for (int k = i+1; k < R; k++)
			{
				T sum = 0;
				for (int j = 0; j < i; j++)
					sum = L[k][j] * U[j][i];
				L[k][i] = (M[k][i] - sum) / U[i][i];
			}
		}
	}

	template <typename T, unsigned int N>
	void lu(const Matrix_t<T, N, N>& M, Matrix_t<T, N, N>& L, Matrix_t<T, N, N>& U)
	{
		for (int i = 0; i < N; i++)
		{
			L[i][i] = 1;

			for (int j = 0; j < i; j++)
				U[i][j] = 0;

			for (int j = i + 1; j < N; j++)
				L[i][j] = 0;

			for (int k = i; k < N; k++)
			{
				T sum = 0;
				for (int j = 0; j < i; j++)
					sum += L[i][j] * U[j][k];
				U[i][k] = M[i][k] - sum;
			}

			for (int k = i + 1; k < N; k++)
			{
				T sum = 0;
				for (int j = 0; j < i; j++)
					sum = L[k][j] * U[j][i];
				L[k][i] = (M[k][i] - sum) / U[i][i];
			}
		}
	}

	template <typename T, unsigned int N>
	void plu(Matrix_t<T, N, N> M, Matrix_t<T, N, N>& P, Matrix_t<T, N, N>& L, Matrix_t<T, N, N>& U)
	{
		Vector_t<int, N> p;
		P = Matrix_t<T, N, N>(0);

		for (int i = 0; i < N; i++)
		{
			p[i] = i;
			L[i][i] = 1;
		}

		for (int k = 0; k < N-1; k++)
		{
			if (unsigned int row = maxloc(abs(column(M, k))); row != k)
				std::swap(p[k], p[row]);

			for (int i = k+1; i < N; i++)
				M[p[i]][k] /= M[p[k]][k];

			for (int i = k+1; i < N; i++)
				for (int j = k+1; j < N; j++)
					M[p[j]][i] -= M[p[j]][k] * M[p[k]][i];

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < i; j++)
				{
					L[i][j] = M[p[i]][j];
					U[i][j] = 0;
				}

				for (int j = i; j < N; j++)
					U[i][j] = M[p[i]][j];

				for (int j = i+1; j < N; j++)
					L[i][j] = 0;
			}
		}

		for (int i = 0; i < N; i++)
			P[p[i]][i] = 1;
	}

	template <typename T, unsigned int N>
	Matrix_t<T, N, N> cholesky(const Matrix_t<T, N, N>& M)
	{
		Matrix_t<T, N, N> res(0);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				T sum = 0;
				for (int k = 0; k < j; k++)
					sum += res[i][k] * res[j][k];
				res[i][j] = (M[i][j] - sum) / res[j][j];
			}

			T sum = 0;
			for (int j = 0; j < i; j++)
				sum += res[i][j] * res[i][j];
			res[i][i] = std::sqrt(M[i][i] - sum);
		}

		return res;
	}
}