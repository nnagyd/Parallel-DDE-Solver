#ifndef PARALLELDDE
#define PARALLELDDE

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "C:/Users/nnagy/Documents/Egyetem/HDS/VCL/version2-2.01.02/vectorclass.h"
#include "ParalellDDE_Common.h"
#define vecSize 4

template<unsigned int nrOfVars, unsigned int nrOfDelays = 0, unsigned int nrOfDenseVars = 0, unsigned int nrOfParameters = 0, unsigned int nrOfEvents = 0, unsigned int nrOfUnroll = 1>
class ParalellDDE
{
private:
	//sizes
	unsigned int nrOfInitialPoints, nrOfSteps, nrOfExtraPoints; //user given
	unsigned int denseOutputMemoryLength, memoryId;//calculated

	//integration memory to every variable with dense output
	double *** __restrict xVals, *** __restrict xdVals, * __restrict tVals;
	/*
	Integration memory layout
	length of outer dimension: denseOutputMemorySize
	length of middle dimension: nrOfDenseVars
	length of inner dimension: nrOfUnroll * vecSize
	index in inner dimension: unroll * vecSize + vecId
	*/

	//mesh points
	unsigned int meshLen; //number of mesh points
	unsigned int meshId; //index of next mesh point
	double* mesh; //synchronized mesh points for variables -> avoid divergence
	double meshPrecision;
	int* meshType; //1 = C1 disc or lower (two points saved)	2 = C2 disc. or higher (one point saved, stepped on that point)

	//delays
	double t0[nrOfDelays];
	unsigned int delayIdToDenseId[nrOfDelays]; //lookup table:	delay -> dense variable index
	unsigned int varIdToDenseId[nrOfVars]; //lookup table: variable -> dense variable index
	unsigned int denseIdToVarId[nrOfDenseVars]; //lookup table: dense variable -> variable index

	//dynamic integration variables
	Vec4d x[nrOfVars * nrOfUnroll];
	double t, dt; //synchronized for all threads

	//initial condition
	Vec4d x0[nrOfVars * nrOfUnroll];

	//temporary integration variables
	Vec4d kAct[nrOfVars * nrOfUnroll], kSum[nrOfVars * nrOfUnroll], xTmp[nrOfVars * nrOfUnroll], xDelay[nrOfDelays * nrOfUnroll];
	double tTmp;
	/*
	integration variable layout
	indexing: unroll*nrOfvars + var	
	*/

	//dense output
	Vec4d xb[nrOfDelays * nrOfUnroll], xn[nrOfDelays * nrOfUnroll], xdb[nrOfDelays * nrOfUnroll], xdn[nrOfDelays * nrOfUnroll];
	double tb[nrOfDelays], deltat[nrOfDelays], pdeltat[nrOfDelays];
	unsigned int lastIndex[nrOfDelays];

	//user given integration range
	double tStart, tEnd, dtBase;

	//parameters
	Vec4d * p;

	//initial discontinouities
	double a;
	double* discC0Init, *discC1Init, *discInit;
	unsigned int nrOfC0, nrOfC1, nrOfDisc;
public:
	//constructor
	ParalellDDE() : nrOfExtraPoints(10), meshPrecision(1e-14), meshLen(0), nrOfC0(0), nrOfC1(0), nrOfDisc(0) 
	{
		for (size_t i = 0; i < nrOfDelays; i++)
		{
			lastIndex[i] = 0;
		}
		p = new Vec4d[nrOfParameters * nrOfUnroll];
	}; 

	//destruktor
	~ParalellDDE() {}; //does nothing

	//set functions
	void setIntegrationRange(double tStart, double tEnd)
	{
		this->tStart = tStart;
		this->tEnd = tEnd;
	}
	void setStepSize(double dt)
	{
		this->dtBase = dt;
		this->nrOfSteps = unsigned((tEnd - tStart) / dt);
	}
	void setNrOfInitialPoints(unsigned int nrOfInitialPoints)
	{
		this->nrOfInitialPoints = nrOfInitialPoints;
	}
	void setNrOfSteps(unsigned int nrOfSteps)
	{
		this->nrOfSteps = nrOfSteps;
		this->dtBase = (tEnd - tStart) / double(nrOfSteps);
	}
	void setNrOfDisc(unsigned int nrOfDisc)
	{
		this->nrOfDisc = nrOfDisc;
	}
	void setNrOfExtraPoints(unsigned int nrOfExtraPoints)
	{
		this->nrOfExtraPoints = nrOfExtraPoints;
	}
	void setInitialTvals(double* tVals)
	{
		//allocate tVals
		this->denseOutputMemoryLength = nrOfInitialPoints + nrOfSteps + meshLen * 2 + nrOfExtraPoints;
		this->tVals = new double[denseOutputMemoryLength];

		//allocate xVals and xdVals
		this->xVals = new double** [denseOutputMemoryLength];
		this->xdVals = new double** [denseOutputMemoryLength];
		for (size_t i = 0; i < denseOutputMemoryLength; i++)
		{
			this->xVals[i] = new double* [nrOfDenseVars];
			this->xdVals[i] = new double* [nrOfDenseVars];
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				this->xVals[i][j] = new double[nrOfUnroll * vecSize];
				this->xdVals[i][j] = new double[nrOfUnroll * vecSize];
			}

		}

		for (size_t i = 0; i < nrOfInitialPoints; i++)
		{
			this->tVals[i] = tVals[i];
		}
	}
	void setInitialXvals(double* xVals, double* xdVals, unsigned int xVarIndex, unsigned int xDenseVarIndex)
	{
		this->varIdToDenseId[xVarIndex] = xDenseVarIndex;
		this->denseIdToVarId[xDenseVarIndex] = xVarIndex;
		for (size_t i = 0; i < nrOfInitialPoints; i++)
		{
			for (size_t j= 0; j < nrOfUnroll * vecSize; j++)
			{
				this->xVals[i][xDenseVarIndex][j] = xVals[i];
				this->xdVals[i][xDenseVarIndex][j] = xdVals[i];
			}
			
		}
	}
	void setInitialDisc(double* C1disc, unsigned int nrOfC1, double* C0disc, unsigned int nrOfC0)
	{
		this->nrOfC0 = nrOfC0;
		this->nrOfC1 = nrOfC1;
		this->nrOfDisc = nrOfC0 + nrOfC1;

		this->discC0Init = new double[nrOfC0];
		this->discC1Init = new double[nrOfC1];
		this->discInit = new double[nrOfDisc];

		for (size_t i = 0; i < nrOfC0; i++)
		{
			this->discC0Init[i] = C0disc[i];
			this->discInit[i] = C0disc[i];
		}

		for (size_t i = 0; i < nrOfC1; i++)
		{
			this->discC1Init[i] = C1disc[i];
			this->discInit[i + nrOfC0] = C1disc[i];
		}

	}
	void setMesh(double* mesh, int* meshType, unsigned int meshLen)
	{
		this->meshLen = meshLen;
		for (size_t i = 0; i < meshLen; i++)
		{
			this->mesh[i] = mesh[i];
			this->meshType[i] = meshType[i];
		}
	}
	void setMeshPrecision(double prec)
	{
		this->meshPrecision = prec;
	} //should be at max 14 order smaller than the order of the step times
	void setRange(double tStart, double tEnd)
	{
		this->tStart = tStart;
		this->tEnd = tEnd;
	}
	void addDelay(unsigned int delayId,double t0, unsigned int denseVarId)
	{
		this->t0[delayId] = t0;
		this->delayIdToDenseId[delayId] = denseVarId;
	}
	void setInitialConditions(Vec4d* x0) {};
	void setX0(double x0)
	{
		for (size_t i = 0; i < nrOfUnroll * nrOfVars; i++)
		{
			this->x0[i] = x0;
		}
	}
	void addMeshPoint(double t, int type)
	{
		unsigned int newMeshLen = meshLen + 1;
		double* newMesh = new double[newMeshLen];
		int* newMeshType = new int[newMeshLen];

		bool placed = false;
		for (size_t i = 0; i < meshLen; i++)
		{
			if (t < mesh[i])
			{
				if (!placed)
				{
					newMesh[i] = t;
					newMeshType[i] = type;
					placed = true;
					i++;
				}
				newMesh[i] = mesh[i-1];
				newMeshType[i] = meshType[i-1];
			}
			else
			{
				newMesh[i] = mesh[i];
				newMeshType[i] = meshType[i];
			}
			
		}

		if (!placed)
		{
			newMesh[meshLen] = t;
			newMeshType[meshLen] = type;
		}
		else
		{
			newMesh[meshLen] = mesh[meshLen - 1];
			newMeshType[meshLen] = meshType[meshLen - 1];;
		}
		delete mesh, meshType;
		meshLen = newMeshLen;
		mesh = newMesh;
		meshType = newMeshType;
	}
	void setParameters(double* p, unsigned int parameterId) //size of p should be nrOfUnroll * nrOfParameter * 4
	{
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			this->p[nrOfParameters * i + parameterId].load(p + 4 * i);
		}
	}
	void setParameters(double p, unsigned int parameterId) //set all parameters to the same
	{
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			this->p[nrOfParameters * i + parameterId] = p;
		}
	}

	//get functions
	unsigned int getDenseOutputMemoryLength() //returns memory length
	{
		return denseOutputMemoryLength;
	}
	unsigned int getDenseOutputMemorySize() //returns whole memory size in bytes
	{
		return 2 * denseOutputMemoryLength * nrOfDenseVars * nrOfUnroll * vecSize * sizeof(double) + denseOutputMemoryLength * sizeof(double);
	}
	double* getMesh()
	{
		double* meshTmp = new double[meshLen];
		for (size_t i = 0; i < meshLen; i++)
		{
			meshTmp[i] = mesh[i];
		}
		return meshTmp;
	}
	double* getMeshType()
	{
		int* meshTmp = new double[meshLen];
		for (size_t i = 0; i < meshLen; i++)
		{
			meshTmp[i] = meshType[i];
		}
		return meshTmp;
	}
	unsigned int getMeshLen()
	{
		return meshLen;
	}
	double* getEndValues(int varId = 0)
	{
		double* endVals = new double[vecSize * nrOfUnroll];

		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			for (size_t j = 0; j < vecSize; j++)
			{
				int id = i * vecSize + j;
				endVals[id] = x[nrOfVars * i + varId][j];
			}
		}
		return endVals;
	}

	//calculate functions - alternative to set, it does the job itself
	void calculateIntegrationMesh()
	{
		double* discC1Full = new double[nrOfC1 + 1];
		discC1Full[0] = 0.0;
		for (size_t i = 1; i < nrOfC1 + 1; i++)
		{
			discC1Full[i] = discC1Init[i - 1];
		}
		calculateFullIntegrationMesh(&mesh, &meshType, &meshLen, discC0Init, nrOfC0, discC1Full, nrOfC1 + 1, t0, nrOfDelays, tStart, meshPrecision);
		this->addMeshPoint(tEnd, 1);
		delete discC1Full;
	}
	void calculateInitialTvals(double t0max)
	{
		double *tTmp = linspaceDisc(t0max, tStart, nrOfInitialPoints, discInit, nrOfDisc);
		this->setInitialTvals(tTmp);
	}
	void calculateInitialXvals(double x0(double t, int dir), double xd0(double t, int dir), unsigned int xVarIndex, unsigned int xDenseVarIndex)
	{
		double* xTmp = discretize(x0, tVals, nrOfInitialPoints);
		double* xdTmp = discretize(xd0, tVals, nrOfInitialPoints);
		this->setInitialXvals(xTmp, xdTmp, xVarIndex, xDenseVarIndex);
	}

	//integration functions
	void integrate(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p));
	void RK4step(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p));

	//dense output
	void calculatedenseOutput(double atT, int delayId);
	void calculateAllDelay(double atT);

	//save functions
	void save(std::string filename, unsigned int varId = 0,int precision = 16) 
	{
		std::ofstream ofs(filename);
		if (!ofs.is_open())
		{
			std::cout << "Cannot open file: " << filename << std::endl;
			return;
		}
		ofs << std::setprecision(precision);
		for (size_t i = nrOfInitialPoints; i <= memoryId; i++)
		{
			ofs << tVals[i] << "\t";
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				ofs << xVals[i][j][varId] << "\t";
			}
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				ofs << xdVals[i][j][varId] << "\t";
			}
			ofs << std::endl;
		}
	}

	//debug functions
	void printMesh(int precision = 6)
	{
		std::cout << std::setprecision(precision) << "Mesh Length = " << meshLen << std::endl;
		for (size_t i = 0; i < meshLen; i++)
		{
			std::cout << mesh[i] << "\tType: " << meshType[i] << std::endl;
		}
	}
	void printInitialPoints(int precision = 6)
	{
		std::cout << std::setprecision(precision) << "Number of initial points = " << nrOfInitialPoints << std::endl;
		for (size_t i = 0; i < nrOfInitialPoints; i++)
		{
			std::cout << i << "\tt = " << std::setw(precision + 6) << tVals[i] << "\t x = ";
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				for (size_t k = 0; k < vecSize * nrOfUnroll; k++)
				{
					std::cout << std::setw(precision + 6) << xVals[i][j][k] << "\t";
				}
			}
			std::cout << "xd = ";
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				for (size_t k = 0; k < vecSize * nrOfUnroll; k++)
				{
					std::cout << std::setw(precision + 6) << xdVals[i][j][k] << "\t";
				}
			}
			std::cout << std::endl;
		}
	}
	void printParameters(int precision = 6)
	{
		std::cout << std::setprecision(precision) << "Number of parameters = " << nrOfParameters << std::endl;

		for (size_t i = 0; i < nrOfParameters; i++)
		{
			std::cout << "Parameter 1:" << std::endl;
			for (size_t j = 0; j < nrOfUnroll; j++)
			{
				for (size_t k = 0; k < vecSize; k++)
				{
					std::cout << p[j * nrOfParameters + i][k] << std::endl;
				}
			}
			std::cout << std::endl;
		}
	}
};

#endif