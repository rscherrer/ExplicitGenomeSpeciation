This folder contains debug and release builds of the program. Those were compiled and linked in QtCreator on Linux.

The Qmake configuration files used are in the root folder, EGS_test.pro (debug) and EGS.pro (release).

The debug version implements tests. The executable is EGS_test.
The release version is optimized for speed. The executable is EGS.

You can move or copy the executables somewhere else on your computer to run them.

If an executable does not run on your platform, you can recompile the program with the build system of your choice from the source code present in ../library (or ../test for the debug version).

# Run the tests

To run the tests (i.e. the debug version), go where the EGS_test executable is, and run:

```{bash}
./EGS_test
```

The debug version should print to the screen the results of its tests as it runs.

# Run a simulation

To run a simulation with default parameters, go where the EGS executable is and run:

```{bash}
./EGS
```

You can have a look at the default parameters in file /parameters_default.txt. Otherwise, they are specified in the manual.

To run a simulation with custom parameter values, run the program with a command-line argument:

```{bash}
./EGS parameters.txt
```

where parameters.txt is a file containing the parameters. 

If you are going to use parameter files, please familiarize yourself with the structure of, for example, /parameters_example.txt.

Parameters are read form the parameter file expecting a specific syntax:

Lines start with the name of the parameter, followed by a space and the value of the parameter.
For parameters with multiple values, add the extra values after the first one, all separated by spaces.
Invalid parameter names will crash the program.
Invalid parameter values (e.g. negative number of genes) will crash the program.
The parameter file parameters_example.txt gives you an overview of all the parameters that can be user-defined.
New parameters do not actually need a new line, a space is enough.

The release version should print to the screen the status of the simulation (started, ended) and whether the population went extinct.
If parameter record is set to one, the program will also output data files, in .dat format, in the working directory the program was run in.
Those .dat files record multiple variables of the simulation through time.
A detailed explanation of the content of these files will come soon.
Those files store floating point numbers,  written in binary.


