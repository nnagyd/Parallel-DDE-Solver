/*
Example 3.: Logistic equation, parameter sweep
1 variable
1 delay
1 parameter
0 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "source/ParallelDDE_RK4.h"

void logistic1(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = x[0] * (p[0] - delay[0]);
}

double x0(double t, int)
{
    return 1.5 - cos(t);
}

double dx0(double t, int)
{
    return sin(t);
}

int main()
{
    //initialize solver
    const unsigned Vars = 1;
    const unsigned Delays = 1;
    const unsigned DenseVars = 1;
    const unsigned Parameters = 1;
    const unsigned Steps = 1000;
    const unsigned InitialPoints = 100;
    const unsigned Unroll = 1;
    ParallelDDE_RK4<Vars, Delays, DenseVars, Parameters, 0, Unroll> solver;

    //parameter sweep
    const unsigned nrOfEquations = 18144;
    double* parameterList = linspace(0, 4, nrOfEquations);

    //basic setup
    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(Steps);
    solver.setNrOfInitialPoints(InitialPoints);

    //delay
    solver.addDelay(0, 1, 0);

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();
    solver.printMesh();

    //initial function
    solver.calculateInitialTvals(-1.0);
    solver.calculateInitialXvals(x0, dx0, 0, 0);

    //initial conditions
    solver.setX0(0.5, 0);

    //save end values
    std::ofstream ofs("example3_endvals.txt");
    ofs << std::setprecision(16);

    //start parameter sweep
    double iStart = seconds();
    for (size_t i = 0; i < nrOfEquations; i += Unroll * vecSize)
    {
        //parameter
        solver.setParameters(parameterList + i, 0);

        //integrate
        solver.integrate(logistic1);

        double* endVals = solver.getEndValues();

        for (size_t j = 0; j < Unroll * vecSize; j++)
        {
            ofs << parameterList[i + j] << "\t" << endVals[j] << std::endl;
        }

        delete endVals;
    }
    double iElapsed = seconds() - iStart;

    //save all values
    solver.save("example3_timeseries.txt");

    //close file
    ofs.flush();
    ofs.close();

    std::cout << "Number of equations= " << nrOfEquations << "\tUnroll= " << Unroll << "\tTime= " << iElapsed << "s\n";
    delete parameterList;
}
