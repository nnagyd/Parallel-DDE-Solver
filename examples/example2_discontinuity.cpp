/*
Example 2.: Basic equation 2
1 variable
1 delay
1 parameter
2 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/
#include <iostream>
#include <string>
#include "source/ParalellDDE_RK4.h"

std::string location = "example2_results.txt";

double x0(double t, int dir)
{
    if (t < -2.0 || (t == -2.0 && dir == -1)) return 0;
    if (t < -1.0 || (t == -1.0 && dir == -1)) return 1;
    return 0;
}

double xd0(double t, int)
{
    return 0;
}


void analitic2(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
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
    const unsigned nrOfSteps = 30;
    ParalellDDE_RK4<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters> solver;

    //basic setup
    solver.setRange(0.0, 25.0);
    solver.setNrOfSteps(nrOfSteps);

    //delay
    solver.addDelay(0, 3, 0);

    //discontinous points
    const unsigned nrOfC0disc = 2;
    double* C0disc = new double[nrOfC0disc]{ -2.0, -1.0 };
    solver.setInitialDisc(NULL, 0, C0disc, nrOfC0disc); //there isn't any C1 discontinous point in the initial condition

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();
    solver.printMesh();

    //parameter
    solver.setParameters(1.0, 0);

    //initial function
    solver.calculateInitialTvals(-3.0);
    solver.calculateInitialXvals(x0, xd0, 0, 0);

    //initial conditions
    solver.setX0(0.0, 0);

    //integrate
    solver.integrate(analitic2);

    //save
    solver.save(location);
}
