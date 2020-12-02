#pragma once

#include <unordered_map>
#include <string>

template <typename T, unsigned int R, unsigned int C>
class Sparse_t
{
public:
	std::unordered_map<unsigned int, T> data[R];

private:
	Sparse_t& add(const Sparse_t& m)
	{
		for (int i = 0; i < R; i++) {
			for (const auto e : m.data[i]) {
				auto& ref = data[i];
				if (const auto it = ref.find(e.first); it != ref.end())
					it->second += e.second;
				else
					ref.insert(e);
			}
		}

		return *this;
	}

	Sparse_t& sub(const Sparse_t& m)
	{
		for (int i = 0; i < R; i++) {
			for (const auto e : m.data[i]) {
				auto& ref = data[i];
				if (const auto it = ref.find(e.first); it != ref.end())
					it->second -= e.second;
				else
					ref[e.first] = -e.second;
			}
		}

		return *this;
	}

	Sparse_t& mul(const Sparse_t& m)
	{
		Sparse_t res;

		for (int i = 0; i < R; i++)
		{
			for (int j = 0; j < C; j++)
			{
				T sum = 0;
				for (auto& e : data[i]) {
					const auto& ref = m.data[e.first];
					if (const auto it = ref.find(j); it != ref.end())
						sum += e.second * it->second;
				}

				res.data[i][j] = sum;
			}
		}

		*this = res;

		return *this;
	}

	Sparse_t& mul(const T& s)
	{
		for (int i = 0; i < R; i++)
			for (auto& e : data[i])
				e.second *= s;

		return *this;
	}

	Sparse_t& div(const T& s)
	{
		for (int i = 0; i < R; i++)
			for (auto& e : data[i])
				e.second /= s;

		return *this;
	}

public:
	Sparse_t() {}

	Sparse_t(T* arr)
	{
		for (int i = 0; i < R; i++)
			for (int j = 0; j < C; j++)
				if (const T t = arr[i * C + j]; t)
					data[i][j] = t;
	}

	Sparse_t(const Sparse_t& m)
	{
		for (int i = 0; i < R; i++)
			data[i] = m.data[i];
	}

	friend Sparse_t operator+(Sparse_t lhs, const Sparse_t& rhs)
	{
		return lhs.add(rhs);
	}

	friend Sparse_t operator-(Sparse_t lhs, const Sparse_t& rhs)
	{
		return lhs.sub(rhs);
	}

	friend Sparse_t operator*(Sparse_t lhs, const Sparse_t& rhs)
	{
		return lhs.mul(rhs);
	}

	friend Sparse_t operator*(Sparse_t m, const T& s)
	{
		return m.mul(s);
	}

	friend Sparse_t operator*(const T& s, Sparse_t m)
	{
		return m.mul(s);
	}

	friend Sparse_t operator/(Sparse_t m, const T& s)
	{
		return m.div(s);
	}

	Sparse_t& operator+=(const Sparse_t& m)
	{
		return add(m);
	}

	Sparse_t& operator-=(const Sparse_t& m)
	{
		return sub(m);
	}

	Sparse_t& operator*=(const Sparse_t& m)
	{
		return mul(m);
	}

	Sparse_t& operator*=(const T& s)
	{
		return mul(s);
	}

	Sparse_t& operator/=(const T& s)
	{
		return div(s);
	}

	bool operator==(const Sparse_t& m) const
	{
		for (int i = 0; i < R; i++)
		{
			const auto& ref1 = data[i];
			const auto& ref2 = m.data[i];

			if (ref1.size() != ref2.size())
				return false;

			for (const auto& e : ref1)
			{
				if (const auto it = ref2.find(e.first); it != ref2.end()) {
					if (e.second != it->second) { return false; }
				} else
					return false;
			}
		}

		return true;
	}

	bool operator!=(const Sparse_t& m) const
	{
		for (int i = 0; i < R; i++)
		{
			const auto& ref1 = data[i];
			const auto& ref2 = m.data[i];

			if (ref1.size() != ref2.size())
				return true;

			for (const auto& e : ref1)
			{
				if (const auto it = ref2.find(e.first); it != ref2.end()) {
					if (e.second != it->second) { return true; }
				}
				else
					return true;
			}
		}

		return false;
	}

	T sparsity() const
	{
		unsigned int non_zero_elems = 0;
		for (int i = 0; i < R; i++)
			non_zero_elems += data[i].size();

		return static_cast<T>(1) - static_cast<T>(non_zero_elems) / (R * C);
	}

	std::string toString() const
	{
		std::string res = "Sparse_t<";
		res += std::string(typeid(T).name()) + ", " + std::to_string(R) + ", " + std::to_string(C) + ">{";

		for (int i = 0; i < R; i++) {
			const auto it = data[i].find(0);
			res += "\n\t[" + ((it != data[i].end()) ? std::to_string(it->second) : std::to_string(static_cast<T>(0)));

			for (int j = 1; j < C; j++) {
				const auto it = data[i].find(j);
				res += ", " + ((it != data[i].end()) ? std::to_string(it->second) : std::to_string(static_cast<T>(0)));
			}

			res += "]";
		}
		res += "\n}";

		return res;
	}

	friend std::ostream& operator<<(std::ostream& os, const Sparse_t<T, R, C>& m)
	{
		os << m.toString();
		return os;
	}
};