/*
Example 1.: Basic equation 1
1 variable
1 delay
1 parameter
0 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/
#include <iostream>
#include <string>
#include "source/ParalellDDE_RK4.h"

std::string location = "example1_results.txt";

double x0(double t, int)
{
    return 1;
}

double xd0(double t, int)
{
    return 0;
}


void analitic1(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = -p[0] * delay[0];
}

int main()
{
    //initialize solver
    const unsigned nrOfVars = 1;
    const unsigned nrOfDelays = 1;
    const unsigned nrOfDenseVars = 1;
    const unsigned nrOfParameters = 1;
    const unsigned nrOfSteps = 100;
    ParalellDDE_RK4<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters> solver;

    //basic setup
    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(nrOfSteps);

    //delay
    solver.addDelay(0, 1, 0);

    //parameter
    solver.setParameters(1.0, 0);

    //initial function
    solver.calculateInitialTvals(-2);
    solver.calculateInitialXvals(x0,xd0,0,0);

    //initial conditions
    solver.setX0(1.0, 0);

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();

    //integrate
    solver.integrate(analitic1);

    //save
    solver.save(location);
}
