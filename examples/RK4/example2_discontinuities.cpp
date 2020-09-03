/*
Example 2.: Basic equation 2
1 variable
1 delay
1 parameter
2 discontinous points

Created by: Daniel Nagy, 24/08/2020
*/

#include "source/ParallelDDE_RK4.h"

void analitic2(Vec4d* xd, double t, Vec4d* x, Vec4d* delay, Vec4d* p)
{
    xd[0] = -p[0] * delay[0];
}

double x0(double t, int dir)
{
    if (t < -2.0 || (t == -2.0 && dir == -1)) return 0;
    if (t < -1.0 || (t == -1.0 && dir == -1)) return 1;
    return 0;
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
    const unsigned Parameters = 1;
    ParallelDDE_RK4<Vars, Delays, DenseVars, Parameters> solver;

    //basic setup
    solver.setRange(0.0, 25.0);
    solver.setNrOfSteps(235);

    //delay
    solver.addDelay(0, 3, 0);

    //discontinous points
    const unsigned nrOfC0disc = 2;
    double* C0disc = new double[nrOfC0disc] { -2.0, -1.0 };
    solver.setInitialDisc(NULL, 0, C0disc, nrOfC0disc); //there isn't any C1 discontinous point in the initial condition

    //mesh
    solver.calculateIntegrationMesh();
    solver.printMesh();

    //parameter
    solver.setParameters(1.0, 0);

    //initial function
    solver.calculateInitialTvals(-3.0);
    solver.calculateInitialXvals(x0, dx0, 0, 0);

    //initial conditions
    solver.setX0(0.0, 0);

    //integrate
    solver.integrate(analitic2);

    //save
    solver.save("example2_results.txt");

    //delete allocated memory
    delete C0disc;
}
