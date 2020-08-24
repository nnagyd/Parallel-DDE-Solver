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
#include "source/ParalellDDE_RK4.h"

double x0(double t, int)
{
    return 1.5 - cos(t);
}

double xd0(double t, int)
{
    return sin(t);
}


void logistic1(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = x[0] * (p[0] - delay[0]);
}

int main()
{
    //initialize solver
    const unsigned nrOfVars = 1;
    const unsigned nrOfDelays = 1;
    const unsigned nrOfDenseVars = 1;
    const unsigned nrOfParameters = 1;
    const unsigned nrOfSteps = 10000;
    const unsigned unroll = 4;
    ParalellDDE_RK4<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters,0,unroll> solver;

    //parameter sweep
    const unsigned parameters = 18144;
    double* parameterList = linspace(0, 4, parameters);

    //basic setup
    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(nrOfSteps);

    //delay
    solver.addDelay(0, 1, 0);

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();

    //initial function
    solver.calculateInitialTvals(-1.0);
    solver.calculateInitialXvals(x0, xd0, 0, 0);

    //initial conditions
    solver.setX0(0.5, 0);

    //save end values
    std::ofstream ofs("endvals.txt");
    ofs << std::setprecision(16);

    //start parameter sweep
    double iStart = seconds();
    for (size_t i = 0; i < parameters; i += unroll * vecSize)
    {
        //parameter
        solver.setParameters(parameterList + i, 0);

        //integrate
        solver.integrate(logistic1);

        double* endVals = solver.getEndValues();

        for (size_t j = 0; j < unroll * vecSize; j++)
        {
            ofs << parameterList[i + j] << "\t" << endVals[j] << std::endl;
        }

        delete endVals;
    }
    double iElapsed = seconds() - iStart;

    //close file
    ofs.flush();
    ofs.close();

    std::cout << "Parameters= " << parameters << "\tUnroll= " << unroll << "\tTime= " << iElapsed << "s\n";
    delete parameterList;
}
