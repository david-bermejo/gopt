#include "math/gopt.hpp"

#include <iostream>
#include <cassert>

#define STR(x) #x
#define ASSERT(x) if (!((x))) { std::cout << "Assertion failed: (" << STR(x) << ")" << std::endl; abort(); }

void unit_tests()
{
	using namespace gopt;

	std::cout << "Vector tests..." << std::endl;

	ASSERT((Vec3f(1, 5, 7) == Vec3f(1, 5, 7)));
	ASSERT((Vec3(1) == Vec3(1, 1, 1)));

	ASSERT((Vec3(1, 4, 5) + Vec3(3) == Vec3(4, 7, 8)));
	ASSERT((Vec3(5, 7, 8) - Vec3(3) == Vec3(2, 4, 5)));
	ASSERT((2 * Vec3(6) == Vec3(12)));
	ASSERT((Vec3(4, 3, 2) * Vec3(5, 6, 4) == 46));
	ASSERT((Vec3(9) / 3 == Vec3(3)));

	ASSERT((cross(Vec3(2, 3, 7), Vec3(7, 3, 1)).almost_equal(Vec3(-18, 47, -15))));

	std::cout << "Quaternion tests..." << std::endl;

	ASSERT(((Quat(pi<>, Vec3(0, 0, 1)) * Vec3(3, 3, 0)).almost_equal(Vec3(-3, -3, 0))));
	ASSERT(((Quat(pi<>/2, Vec3(0, 0, 1)) * Quat(pi<>/2, Vec3(0, 0, 1))).almost_equal(Quat(pi<>, Vec3(0, 0, 1)))));
	ASSERT(((Quaternion_t<float>(pi<float>/2, Vec3f(0, 0, 1)) * Quaternion_t<float>(pi<float>/2, Vec3f(0, 0, 1))).almost_equal(Quaternion_t<float>(pi<float>, Vec3f(0, 0, 1)))));
	ASSERT(((Quat(pi<>/2, Vec3(1, 0, 0)) * Quat(pi<>/2, Vec3(1, 0, 0))).almost_equal(Quat(pi<>, Vec3(1, 0, 0)))));
	ASSERT((rotation(Quat(pi<>/4, Vec3(0, 0, 1)) * Quat(pi<>/6, Vec3(0, 1, 0))).almost_equal(rotation(pi<>/4, Vec3(0, 0, 1)) * rotation(pi<>/6, Vec3(0, 1, 0)))));
	ASSERT((normalize(Quat(2, 2, 2, 2)).almost_equal(Quat(1, 1, 1, 1) / std::sqrt(4))));
	ASSERT((Quat(1, 2, 3, 4).z == 4));
	ASSERT((Quat(1, 2, 3, 4)[3] == 4));
	ASSERT((conjugate(Quat(1, 1, 1, 1)) == Quat(1, -1, -1, -1)));
	ASSERT(((Quat(3, 5, 6, 7) * inverse(Quat(3, 5, 6, 7))).almost_equal(Quat(1, 0, 0, 0))));

	std::cout << "All tests passed!" << std::endl;
}