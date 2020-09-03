/*
Example 1.: Basic equation 1
1 variable
1 delay
1 parameter
0 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/
#include <iostream>
#include "source/ParallelDDE_RK4.h"

void analitic1(Vec4d* dx, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    dx[0] = -p[0] * delay[0];
}


double x0(double t, int)
{
    return 1;
}

double dx0(double t, int)
{
    return 0;
}

int main()
{
    //initialize solver
    const unsigned Vars = 1;
    const unsigned Delays = 1;
    const unsigned DenseVars = 1;
    const unsigned ParametersPerEquation = 1;
    ParallelDDE_RK4<Vars, Delays,DenseVars,ParametersPerEquation> solver;

    //basic setup
    solver.setRange(0.0, 10.0);
    solver.setNrOfSteps(100);

    //delay
    solver.addDelay(0, 1, 0);

    //initial function
    solver.calculateInitialTvals(-1);
    solver.calculateInitialXvals(x0, dx0, 0, 0);

    //mesh (only discontinous point is t = 0)
    solver.calculateIntegrationMesh();
    solver.printMesh();

    //initial conditions
    solver.setX0(1.0, 0);

    //parameter
    solver.setParameters(1.0, 0);

    //integrate
    solver.integrate(analitic1);

    //save
    solver.save("example1_results.txt");
}
