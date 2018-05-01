# runge-kutta-4

Explicit 4th order Runge-Kutta solver for scalar and vector ODEs, C and C++ compliant.

## Getting Started

The solver is a couple of C/C++ compliant files, one source and one header. As a result, the source file (<name>.cpp) can be renamed (<name>.c) freely, and used in both languages. The solver is therefore distributed as raw source code that needs to be compiled.
  
The key point is in the following functions:

```
void rungekutta4()
void rungekutta4_vec()
```
As they are the "actual solver", with the closest syntax as MATLAB odexx. Some extras have been added as the integration step h, and the dimensionality of the problem.

The solver yields a txt file with vector or matrix output, as it can be seen in the /tests/ folder.

### Prerequisites

To implement the solver, a standard C/C++ IDE and compiler needs to be installed on the operative system.  To name a few, Visual Studio 20xx (Visual C++) for Windows or Emacs (GCC) for Linux.

### Installing

1. Follow the instructions for downloading the IDE (Visual Studio, Emacs *et similia*) of your preference, and then build a new C or C++ project. 

2. In such project, simply copy and paste the source files from /source/. Remember to paste the rk4solver.c or kr4solver.cpp file in the "resources" folder and all the rk4solver.h files in the "header" folder of the current project.

3. Compile your project either in debug or release, your choice. The rk4main.c file is an example implementation of the solver. The solver's function set is called by:
```
 #include"rk4solver.h"
```
and the rest is just declaration and definition of objective functions and solver parameters. The latter part is no different from any common MATLAB file when using odexx.

Once the code is compiled, you are good to go!

## Running the tests

*No need to think of new functions for testing!*

Two test functions are already provided within the main cpp file, in the form of standard C/C++ functions. Such functions are also nonlinear, aiming to be worst-case scenarios of nonstiff systems.

## Deployment

The solver can be deployed within any C/C++ IDE and compiler, and virtually to any machine that support the former tools. Data analysis is performed by any tool that can process a text file, for instance Gnuplot. 

For this instance, I added a MATLAB file in /source/ to exemplify the data retreival and validation procedure.

## Contributing and Improvements

*Contributing and redistributing the code is free!*

You may choose to modify and improve it to your liking. Possible improvements may include:

1. Increase the accuracy. Can the integration step be modified at runtime, or be adaptive?

2. Ease the computational load. There are many malloc calls, is it a way to reduce them?

3. Expand to stiff systems. Does the Van der Pol oscillator test work here?

But again, they are limited only by imagination.

## Authors

* **Andrea Giuseppe Landella** - [alandella](https://github.com/alandella)

## License

This project is licensed under The GNU General Public License v3.0 - see the [LICENSE.md](https://github.com/alandella/runge-kutta-4/blob/master/LICENSE) file for details.

## Acknowledgments

* Inspiration
* Curiosity
* Boredom

The latter may imply the former two sometimes. Or vice versa, idk.
