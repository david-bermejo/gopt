GOPT
====

GOPT is a C++ template library for defining and solving calculus and linear algebra problems, as well as system dynamics involving differential equations.

Motivation
----------

This project was inspired by Fortran ease with handling vector and matrices operations. Its aim is to provide tools with an easy and clear syntax to solve any linear and non-linear system of differential equations, proper of dynamical systems.  
It is designed being lightweight and speed its main characteristics, so that it can be executed by a microcontroller in a reasonable amount of time.  
The final goal is to provide mathematical tools required to perform guidance and trajectory control of a real dynamical system.

Installation
------------

GOPT is a header-only library, so download it next to the project you want to use with, and it will be ready.  
This project is ensured to be compiled using the MSVC Windows toolchain. Please, communicate if it can be built under any other environment.

Usage
-----

A brief introduction of this library is provided. It is divided into different sections. [vector.hpp](math/vector.hpp), [matrix.hpp](math/matrix.hpp) and [quaternion.hpp](math/quaternion.hpp) give an implementation of the respective objects, while [algorithms.hpp](math/algorithms.hpp) collects functions and algorithms applied to those mathematical objects.  
Including a single file [gopt.hpp](math/gopt.hpp) gives full access to the library.
  
#### 1. Vectors
```c++
#include "math/gopt.hpp"
using namespace gopt;

int main()
{
    Vec3 u(1, 2, 3); // Same as Vector_t<double,3>.
    Vec3 v(2, 4, 1);
    std::cout << cross(u,v) << std::endl;

    return 0;
}
```
```console
Vector<double,3>{-10, 5, 0}
```
  
#### 2. Matrices
```c++
#include "math/gopt.hpp"
using namespace gopt;

int main()
{
    Mat3 A(1, 2, 5, 6, 8, 1, 3, 3, 2); // Same as Matrix_t<double,3,3>.
    Vec3 b(1, 5, 4);
    std::cout << A*b << std::endl;
    std::cout << "Element i=3, j=3: " A[2][2] << std::endl; // Access the element at 3rd row and 3rd column.

    return 0;
}
```
```console
Vector<double,3>{31, 50, 26}
Element i=3, j=3: 2
```
  
#### 3. Quaternions
Unit [quaternions](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) are mathematical objects used for representing space orientations and rotations of objects in three dimensions.  
GOPT also provides a simple yet useful implementation of this object. It includes a set of functions to describe body rotations as a replacement of the matrix form, which is known to suffer from [gimbal lock](https://en.wikipedia.org/wiki/Gimbal_lock).  
  
```c++
#include "math/gopt.hpp"
using namespace gopt;

int main()
{
    Quat q0(1); // Unit quaternion, same as Quat(1,0,0,0).
    std::cout << q0 << std::endl;

    Quat q(pi<>, Vec3(0, 0, 1)); // Define a rotation of 180 degrees around the z axis.
    Vec3 p(1,1,0);
    
    // Note: rotate(Quat, Vec3) implements fast vector rotation algorithm.
    std::cout << rotate(q, p) << std::endl;

    return 0;
}
```
```console
Quaternion<double>{1, 0, 0, 0}
Vector<double,3>{-1, -1, 0}
```
  
#### 4. ODE Solvers
The library also includes a set of [solvers](math/solvers.hpp) for differential equations. They are capable of solving n-dimensional systems of differential equations. At the moment, two of them are implemented:
* 5th order Runge-Kutta-Nystrom (RKN5)
* 4th order Dormand-Prince (DP45)
  
The second one should be preferred, it is faster and less prone to errors. There are substantial differences between them. The first one is used to solve genera 2nd order differential equations in the form: y" = f(y, y', t), while Dormand-Prince algorithm solves first order differential equations in the form: y'=f(y,t).  
For solving a 2nd order differential equation using the latter algorithm, it should be rewritten in the following form:
```console
dy/dt = y'
d2y/dt2 = f(y, y', t)
```
  
In the following section, an example of how to solve the Van der Pol equation is provided.

```c++
#include "math/gopt.hpp"
using namespace gopt;

/**************************************************
 * Test Dormand-Prince method (DP45) to solve
 * the 2nd order Van der Pol differential equation:
 * y" + y' = 0; y(0) = 1; y'(0) = 1
 *
 * Analytical solution:
 * y(t) = sin(t) + cos(t)
 **************************************************/

Vec2 f(const Vec2& x, const double t)
{
    return Vec2{ x[1], -x[0] };
}

double van_der_pol(const double t)
{
    return std::sin(t) + std::cos(t);
}

int main()
{
    const Vec2 x0(1,1);
    const double t0 = 0;
    const double tf = 5;

    // Solve the differential equation
    Vec2 result_def = DP45(f, x0, t0, tf); // Result with default tolerances
    Vec2 result = DP45(f, x0, t0, tf, 1e-9, 1e-9); // Result with more strict tolerances

    std::cout << "Results at t=5:" << std::endl;
    std::cout << "Numerical (default tolerance): " << result_def[0] << std::endl;
    std::cout << "Numerical (stricter tolerance): " << result[0] << std::endl;
    std::cout << "Analytical: " << van_der_pol(tf) << std::endl;

    return 0;
}
```
```console
Results at t=5:
Numerical (default tolerance): -0.674519
Numerical (stricter tolerance): -0.675262
Analytical: -0.675262
```
  
#### 4. Additional examples

[main.cpp](main.cpp) gives an implementation of a constant mass model rocket with trajectory control via thrust vectoring. Two servo motors are used to deflect the main engine to provide yaw and pitch control.  
An additional [Extended Kalman Filter](https://en.wikipedia.org/wiki/Extended_Kalman_filter) ([ekf.hpp](math/ekf.hpp)) is implemented to provide optimal control with the aid of an [Inertial Measurement Unit](https://en.wikipedia.org/wiki/Inertial_measurement_unit), which gives updated readings of the acceleration and angular velocity of the rocket, to get an estimation of the position and attitude of it.  
A PID controller or similar may be used to allow trajectory control.
  
Contributions
-------------

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

License
-------
[MIT](LICENSE.txt)