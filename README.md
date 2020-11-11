# numerical-methods

This project is just an example of how numerical methods were implemented to solve chemical engineering problems in one of my courses.

The files in this project:

"Project 4 Guidelines.pdf" - This file contains the project guidelines for the course. It specifies a continuous stirred tank reactor (CSTR) with some set parameters. The project is intended to determine the temperature of the exit stream as a funtion of time with an added bonus of determining the temperature when the inlet flow rate is changed with varying valve constants.

"Proj4_2.f90" - This is the FORTRAN program that calculates the exit temperature based on the parameters provided in the problem statement. The program uses both numerical integration and the Runge-Kutta method for calculating the temperature and outputs the data into both a .dat file and uses gnuplot to graph the results. This requires that all the files in the output folders be in the same directory as this program as it uses the files: IntMethod.txt, IntMethod_2.txt, RKMethod.txt, RKMethod_2.txt, IntRKMethod.txt, IntRKMethod_2.txt, and so on, as well as the output .dat files to generate the graphs.

"IntMethod" - This folder contains the output data, graphs, and input files for values calculated by the integration method for part one of the project, to run the program these need to be in the same folder as Proj4_2.f90

"RKMethod" - This folder contains the output data, graphs, and input files for values calculated by the Runge-Kutta method for part one of the project, to run the program these need to be in the same folder as Proj4_2.f90

"IntRKMethod" - This folder contains the output data, graphs, and input files for values calculated by both methods when accounting for the second part of the project (accounting for the change in the inlet flow rate), to run the program these need to be in the same folder as Proj4_2.f90