#ifndef PARALLELDDE_COMMON
#define PARALLELDDE_COMMON

#include <iostream>
#include <chrono>
#include <cmath>

inline void sort(double* lst, unsigned int len)
{
	//bubblesort
	for (size_t i = 1; i < len; i++)
	{
		for (size_t j = 1; j < len; j++)
		{
			if (lst[j] < lst[j - 1]) //swap
			{
				double tmp = lst[j];
				lst[j] = lst[j - 1];
				lst[j - 1] = tmp;
			}
		}
	}
}

inline double* discretize(double f(double t, int dir), double* tVals, unsigned int nrOfPoints)
{
	double* lst = new double[nrOfPoints];

	for (size_t i = 0; i < nrOfPoints; i++)
	{
		int direction = -1;
		if (i >= 1 && tVals[i] == tVals[i - 1]) //change direction by double points
		{
			direction = 1;
		}
		lst[i] = f(tVals[i], direction);
	}
	return lst;
}

inline double* linspace(double a, double b, unsigned int steps)
{
	double* lst = new double[steps];
	double dt = (b - a) / (steps - 1);
	double t = a;
	for (size_t i = 0; i < steps; i++)
	{
		lst[i] = t;
		t += dt;
	}
	return lst;
}

inline double* linspaceDisc(double a, double b, unsigned int steps, double* tDisc = NULL, unsigned int nrOfDisc = 0, double eps = 0.0)
{
	double* lst = new double[steps];
	int* discMask = new int[nrOfDisc];
	//set all element to 0, set to 1 if the i. discontinouity is included 
	for (size_t i = 0; i < nrOfDisc; i++)
	{
		discMask[i] = 0;
	}

	double dt = (b - a) / (steps - 2 * nrOfDisc -1);
	double t = a;
	for (size_t i = 0; i < steps; i++)
	{
		bool set = true;
		for (size_t j = 0; j < nrOfDisc; j++)
		{

			if (!discMask[j] && fabs(tDisc[j] - t) < eps) //discontinuity happens at a point //this function can cause problems, CORRECT 
			{
				lst[i] = tDisc[j] - eps;
				lst[i + 1] = tDisc[j];
				lst[i + 2] = tDisc[j] + eps;
				set = false;
				discMask[j] = 1;
				i += 2;
			}
			else if (!discMask[j] && tDisc[j] < t) //discontinuity far from point
			{
				lst[i] = tDisc[j] - eps;
				lst[i + 1] = tDisc[j] + eps;
				discMask[j] = 1;
				i += 2;
			}
		}
		if (set) lst[i] = t;
		t += dt;
	}
	lst[steps - 1] = b; //last step should be exactly the end of the interval
	return lst;
}

//removes every element from original which is bigger than min, and included in toRemove, with a given tolerance
inline double* filter(double* original, unsigned int* nr, double min, double* toRemove = NULL, unsigned int nrToRemove = 0, double tol = 0.0)
{
	unsigned int count = *nr;
	double* unique = new double[count];
	unsigned int uniqueNr = 0;
	for (size_t i = 0; i < count; i++)
	{
		bool set = true;
		for (size_t j = 0; j < uniqueNr; j++)
		{
			if (abs(original[i] - unique[j]) <= tol) //already in new list
			{
				set = false;
			}
		}

		for (size_t j = 0; j < nrToRemove; j++)
		{
			if (abs(original[i] - toRemove[j]) <= tol) //already in concurrent list
			{
				set = false;
			}
		}

		if (set && original[i] > min) //original[i] not in new list yet
		{
			unique[uniqueNr] = original[i];
			uniqueNr++;
		}
	}

	double* filtered = new double[uniqueNr];
	for (size_t i = 0; i < uniqueNr; i++)
	{
		filtered[i] = unique[i];
	}
	*nr = uniqueNr;
	return filtered;
}

inline double* calculateMesh(double* disc, double* delay, unsigned int sizeDisc, unsigned int sizeDelay, unsigned int recursionDepth = 4, unsigned int dummy = 0)
{
	unsigned int newSize = (sizeDisc - dummy) * sizeDelay + sizeDisc;
	double* newA = new double[newSize];
	for (size_t i = 0; i < sizeDisc; i++)
	{
		newA[i] = disc[i];
	}
	for (size_t i = dummy; i < sizeDisc; i++)
	{
		for (size_t j = 0; j < sizeDelay; j++)
		{
			int id = sizeDisc + (i-dummy)*sizeDelay + j;
			newA[id] = disc[i] + delay[j];
		}
	}
	delete disc;
	if (recursionDepth == 1 || recursionDepth == 0) return newA;
	return calculateMesh(newA, delay, newSize, sizeDelay, recursionDepth - 1,sizeDisc);
}

inline unsigned int powInt(unsigned int a, unsigned int x)
{
	if (x == 0) return 0;
	if (x == 1) return a;
	return a * powInt(a, x - 1);
}

inline unsigned int calculateMeshSize(unsigned int sizeA, unsigned int sizeB, unsigned int depth)
{
	if (sizeB == 1) return sizeA * (depth+1);
	return sizeA * (powInt(sizeB, depth +1) - 1) / (sizeB - 1);
}

//calculates the mesh for the integrator given all discontinous points and all delays
inline void calculateFullIntegrationMesh(double ** integrationMesh, int ** integrationMeshType,unsigned int * integrationMeshLen, double* discC0, unsigned int nrOfC0, double* discC1, unsigned int nrOfC1, double * delays, unsigned int nrOfDelays, double tStart = 0.0,  double tol = 1e-10)
{
	//step 1: find every discontinous point after integration tstart point
	double* doubleStepMesh = calculateMesh(discC0, delays, nrOfC0, nrOfDelays, 1);
	unsigned int nrOfDoubleStepMesh = calculateMeshSize(nrOfC0, nrOfDelays, 1);
	doubleStepMesh = filter(doubleStepMesh, &nrOfDoubleStepMesh, tStart,NULL,0,tol);

	//step 2: concatenate discC1 and doubleStepMesh
	unsigned int totalNrOfDiscC1 = nrOfDoubleStepMesh + nrOfC1;
	double* totalDisC1 = new double[totalNrOfDiscC1];
	for (size_t i = 0; i < totalNrOfDiscC1; i++)
	{
		if (i < nrOfDoubleStepMesh)
		{
			totalDisC1[i] = doubleStepMesh[i];
		}
		else
		{
			totalDisC1[i] = discC1[i - nrOfDoubleStepMesh];
		}
	}

	//step 3: calculate C2, C3, C4 discontinous points
	double* mesh = calculateMesh(totalDisC1, delays, totalNrOfDiscC1, nrOfDelays, 4);
	unsigned int meshLen = calculateMeshSize(totalNrOfDiscC1, nrOfDelays, 4);
	mesh = filter(mesh, &meshLen, tStart, NULL, 0, tol);
	sort(mesh, meshLen);
	
	//step 4: create mesh types
	int* meshType = new int[meshLen];
	for (size_t i = 0; i < meshLen; i++)
	{
		//assuming not in double mesh
		meshType[i] = 1;

		for (size_t j = 0; j < nrOfDoubleStepMesh; j++)
		{
			if (fabs(mesh[i] - doubleStepMesh[j]) < tol)
			{
				//if in double mesh
				meshType[i] = 2;
			}
		}
	}

	//step 5: return values through pointers
	*integrationMesh = mesh;
	*integrationMeshLen = meshLen;
	*integrationMeshType = meshType;
	delete doubleStepMesh;
}

//------------------------------ Timer stuff -------------------------------------------
uint64_t micros()
{
	uint64_t us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::
		now().time_since_epoch()).count();
	return us;
}

double seconds()
{
	return double(micros()) / 1000. / 1000.;
}


#endif // !PARALLELDDE_COMMON