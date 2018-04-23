# runge-kutta-4

Explicit 4th order Runge-Kutta solver for scalar and vector ODEs, C and C++ compliant.

## Getting Started

The solver is a couple of C/C++ compliant files, one source and one header. As a result, the source file (<name>.cpp) can be renamed (<name>.c) freely, and used in both languages. The solver is therefore distributed as raw source code that needs to be compiled.
  
The key point is in the following functions:

```
void rungekutta4()
void* rungekutta4_vec()
```
As they are the "actual solver", almost the same as MATLAB, as the initial and final time, initial condition, and the function handle have to be provided. Some extras have been added as the separation h, or the integration step, and the dimensionality of the problem.

The solver output is a txt file with vector or matrix solutions.

### Prerequisites

To implement the solver, a standard C/C++ IDE and compiler needs to be installed on the operative system.  To name a few, Visual Studio 20xx (Visual C++) for Windows or Emacs (GCC) for Linux.

### Installing

1. Follow the instructions for downloading the IDE (Visual Studio, Emacs *et similia*) of your preference, and then build a new C or C++ project. 

2. In such project, simply copy and paste the source files from my /source/.

3. Compile your project either in debug or release, your choice.

Once the code is compiled, you are good to go!

## Running the tests

*No need to think of new functions for testing!*

Two test functions are already provided within the code, in the form of standard C/C++ functions, and the main cpp file contains applicative examples of the solver. Such functions are also nonlinear, aiming to be worst-case scenarios of nonstiff systems.

## Deployment

The solver can be deployed within any C/C++ IDE and compiler, and virtually to any machine that support the former tools. Data analysis is performed by any tool that can process a text file, for instance Gnuplot. 

## Contributing

*Contributing and redistributing the code is free!*

You may choose to modify and improve it to your liking.

## Authors

* **Andrea Giuseppe Landella** - [alandella](https://github.com/alandella)

## License

This project is licensed under The GNU General Public License v3.0 - see the [LICENSE.md](https://github.com/alandella/runge-kutta-4/blob/master/LICENSE) file for details.

## Acknowledgments

* Inspiration
* Curiosity
* Boredom

The latter implies the former two sometimes.
