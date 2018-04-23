# runge-kutta-4

Explicit 4th order Runge-Kutta solver for scalar and vector ODEs, C and C++ compliant.

## Getting Started

This solver is implemented in C++ and it is backward compatible in C, meaning that it compiles in any of the two languages interpreted by the compiler. As such, it is distributed as raw source code that needs to be compiled.

The solver inputs are almost the same as MATLAB, as the initial and final time, initial condition, and the function handle. Some extras have been added as, the separation h or the integration step and the dimensionality of the problem.

The solver output is a *txt file* with vector or matrix solutions.

### Prerequisites

To implement the solver, a standard C/C++ IDE and compiler needs to be installed on the operative system.  To name a few, Visual Studio 20xx (Visual C++) for Windows or Emacs (GCC) for Linux.

### Installing

TODO

```
TODO
```

## Running the tests

No need to add new functions for testing! 

They are already provided within the code in the form of standard C/C++ functions, and the main contains applicative examples of the solver. Those two functions are ad-hoc nonlinear test functions aiming to be worst-case scenarios of nonstiff systems.

## Deployment

The solver can be deployed within any C/C++ IDE and compiler, and virtually to any machine that support the former tools. Data analysis is performed by any tool that can process a text file, for instance Gnuplot. 

## Contributing

Contributing and redistributing the code is free!

## Authors

* **Andrea Giuseppe Landella** - *Exam work* - [alandella](https://github.com/alandella)

## License

This project is licensed under The GNU General Public License v3.0 - see the [LICENSE.md](https://github.com/alandella/runge-kutta-4/blob/master/LICENSE.md) file for details.

## Acknowledgments

* Mostly Inspiration and Boredom, the latter implies the former sometimes
