/*
Example 4.: Logistic 2. Parameter sweep with discontinouities
1 variable
2 delay
1 parameter
4 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "source/ParallelDDE_RK4.h"

#define PI 3.1415926535897932385

double x0(double t, int dir)
{
    if (t < -1.5 || (t == -1.5 && dir == -1))                           return std::cos(4 * PI * t);
    if (t < -std::sqrt(2.0) || (t == -std::sqrt(2.0) && dir == -1))     return t * t;
    if (t < -1.101 || (t == -1.101 && dir == -1))                       return std::exp(t);
    if (t < -0.5 || (t == -0.5 && dir == -1))                           return 0;
    return t + 0.5;
}

double xd0(double t, int dir)
{
    if (t < -1.5 || (t == -1.5 && dir == -1))                           return -4 * PI * std::sin(4 * PI * t);
    if (t < -std::sqrt(2.0) || (t == -std::sqrt(2.0) && dir == -1))     return 2 * t;
    if (t < -1.101 || (t == -1.101 && dir == -1))                       return std::exp(t);
    if (t < -0.5 || (t == -0.5 && dir == -1))                           return 0;
    return 1.0;
}


void logistic2(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = x[0] * delay[1] * (p[0] - delay[0]);
}

int main()
{
    //initialize solver
    const unsigned nrOfVars = 1;
    const unsigned nrOfDelays = 2;
    const unsigned nrOfDenseVars = 1;
    const unsigned nrOfParameters = 1;
    const unsigned nrOfSteps = 10000;
    const unsigned unroll = 8;
    ParallelDDE_RK4<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters, 0, unroll> solver;

    //parameter sweep
    const unsigned parameters = 18144;
    double* parameterList = linspace(0, 2, parameters);

    //basic setup
    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(nrOfSteps);

    //discontinouities
    const unsigned C0 = 3;
    double* C0disc = new double[C0] {-1.5, -sqrt(2.0), -1.101};
    const unsigned C1 = 1;
    double* C1disc = new double[C1] {-0.5};
    solver.setInitialDisc(C1disc, C1, C0disc, C0);

    //delay
    solver.addDelay(0, 1, 0);
    solver.addDelay(1, 2, 0);

    //mesh (only discontinous point is t = 0)
    solver.setMeshPrecision(1e-13);
    solver.calculateIntegrationMesh();


    //initial function
    solver.setNrOfInitialPoints(200);
    solver.calculateInitialTvals(-2.0);
    solver.calculateInitialXvals(x0, xd0, 0, 0);

    //initial conditions
    solver.setX0(0.5, 0);

    //save end values
    std::ofstream ofs("example4_endvals.txt");
    ofs << std::setprecision(16);

    //start parameter sweep
    double iStart = seconds();
    for (size_t i = 0; i < parameters; i += unroll * vecSize)
    {
        //parameter
        solver.setParameters(parameterList + i, 0);

        //integrate
        solver.integrate(logistic2);

        double* endVals = solver.getEndValues();

        for (size_t j = 0; j < unroll * vecSize; j++)
        {
            ofs << parameterList[i + j] << "\t" << endVals[j] << std::endl;
        }
    }
    double iElapsed = seconds() - iStart;

    //close file
    ofs.flush();
    ofs.close();

    std::cout << "Parameters= " << parameters << "\tUnroll= " << unroll << "\tTime= " << iElapsed << "s\n";
}
