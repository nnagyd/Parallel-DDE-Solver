#ifndef PARALLELDDE
#define PARALLELDDE

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "C:/Users/nnagy/Documents/Egyetem/HDS/VCL/version2-2.01.02/vectorclass.h"
#include "ParalellDDE_Common.cpp"
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
	length of outer dimension:	denseOutputMemorySize
	length of middle dimension: nrOfDenseVars
	length of inner dimension:	nrOfUnroll * vecSize
	index in inner dimension:	unroll * vecSize + vecId
	*/

	//mesh points
	unsigned int meshLen; //number of mesh points
	unsigned int meshId; //index of next mesh point
	double* mesh; //synchronized mesh points for variables -> avoid divergence
	double meshPrecision;
	int* meshType; //1 = C1 disc or lower (two points saved)	2 = C2 disc. or higher (one point saved, stepped on that point)

	//delays
	double * t0;
	unsigned int delayIdToDenseId[nrOfDelays]; //	lookup table:	delay index		-> dense variable index
	unsigned int varIdToDenseId[nrOfVars]; //		lookup table:	variable index	-> dense variable index
	unsigned int denseIdToVarId[nrOfDenseVars]; //	lookup table:	dense variable	-> variable index

	//dynamic integration variables
	Vec4d * __restrict x;
	double t, dt; //synchronized for all threads

	//initial condition
	Vec4d *x0;

	//temporary integration variables
	Vec4d * __restrict kAct, * __restrict kSum, * __restrict xTmp, * __restrict xDelay;
	double tTmp;
	/*
	integration variable layout
	indexing: unroll*nrOfvars + var
	*/

	//dense output
	Vec4d * __restrict xb, * __restrict xn, * __restrict xdb, * __restrict xdn;
	double * __restrict tb, * __restrict deltat, * __restrict pdeltat;
	unsigned int lastIndex[nrOfDelays];

	//user given integration range
	double tStart, tEnd, dtBase;

	//parameters
	Vec4d * __restrict p;

	//initial discontinouities
	double a;
	double* discC0Init, *discC1Init, *discInit;
	unsigned int nrOfC0, nrOfC1, nrOfDisc;

	//event stuff
	Vec4d* prevVals;
	Vec4d* newVals;
	double eventPrecision;
	void (*eventLocation)(Vec4d* lst, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p);
	void (*eventIntervention)(int eventId, int eventDir, Vec4db mask, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p);

public:
	//constructor
	ParalellDDE() : nrOfExtraPoints(10), meshPrecision(1e-14), meshLen(0), nrOfC0(0), nrOfC1(0), nrOfDisc(0), nrOfInitialPoints(50), nrOfSteps(100), dtBase(0.01), tStart(0.0), tEnd(10.0), eventPrecision(1e-10)
	{
		for (size_t i = 0; i < nrOfDelays; i++)
		{
			lastIndex[i] = 0;
		}
		prevVals = new Vec4d[nrOfEvents * nrOfUnroll];
		newVals = new Vec4d[nrOfEvents * nrOfUnroll];
		p = new Vec4d[nrOfParameters * nrOfUnroll];
		x = new Vec4d[nrOfVars * nrOfUnroll];
		x0 = new Vec4d[nrOfVars * nrOfUnroll];
		kAct = new Vec4d[nrOfVars * nrOfUnroll];
		kSum = new Vec4d[nrOfVars * nrOfUnroll];
		xTmp = new Vec4d[nrOfVars * nrOfUnroll];
		xDelay = new Vec4d[nrOfDelays * nrOfUnroll];
		xb = new Vec4d[nrOfDelays * nrOfUnroll];
		xn = new Vec4d[nrOfDelays * nrOfUnroll];
		xdb = new Vec4d[nrOfDelays * nrOfUnroll];
		xdn = new Vec4d[nrOfDelays * nrOfUnroll];
		tb = new double[nrOfDelays];
		deltat = new double[nrOfDelays];
		pdeltat = new double[nrOfDelays];
		t0 = new double[nrOfDelays];
	};

	//destruktor
	~ParalellDDE()
	{
		/*delete p, x, x0, kAct, kSum, xTmp, xDelay;
		delete xb, xn, xdb, xdn, tb, deltat, pdeltat, t0;*/
	}; //does nothing

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
	void setX0(double x0)
	{
		for (size_t i = 0; i < nrOfUnroll * nrOfVars; i++)
		{
			this->x0[i] = x0;
		}
	}
	void setX0(double* x0, unsigned int varId) //size of x0 should be at least nrOfUnroll * nrOfParameter * 4
	{
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			this->x0[nrOfVars * i + varId].load(x0 + 4 * i);
		}
	}
	void setX0(double x0, unsigned int varId) //size of x0 should be at least nrOfUnroll * nrOfParameter * 4
	{
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			this->x0[nrOfVars * i + varId] = x0;
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
	void setParameters(double* p, unsigned int parameterId) //size of p should be at least nrOfUnroll * nrOfParameter * 4
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
	void setT(double t)
	{
		this->t = t;
	}
	void setEventPrecision(double prec)
	{
		this->eventPrecision = prec;
	}
	//event functions
	void addEventLocationFunction(void (*eventLocation)(Vec4d* lst, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		this->eventLocation = eventLocation;
	}
	void addEventInterventionFunction(void (*eventIntervention)(int eventId, int eventDir, Vec4db mask, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		this->eventIntervention = eventIntervention;
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
	double getT()
	{
		return t;
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
	void integrate(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		//initialize variables
		t = tStart; //tStart -> t
		memoryId = nrOfInitialPoints;
		tVals[memoryId] = t;
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			//x0 -> x
			for (size_t j = 0; j < nrOfVars; j++)
			{
				x[i * nrOfVars + j] = x0[i * nrOfVars + j];
			}

			//save dense output points
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				unsigned int varId = denseIdToVarId[j];
				x[i * nrOfVars + varId].store(xVals[memoryId][j] + i * vecSize);
			}
		}
		for (size_t i = 0; i < nrOfDelays; i++)
		{
			lastIndex[i] = 0;
		}
		for (size_t i = 0; i < nrOfEvents; i++)
		{
			prevVals[i] = NAN;
			newVals[i] = NAN;
		}
		meshId = 0;

		//integration
		int stepType = 0;
		//0 normal step
		//1 simple mesh point
		//2 double mesh point
		//3 stepsize decrease because event, doesnt accept step
		//4 acceptable step near event
		//5 event is found, we are after that
		bool nearEvent = false;
		while (t < tEnd)
		{
			if (stepType == 0) //if no double point neccessary
			{
				dt = dtBase;
			}
			if (stepType == 1) //if simple mesh step was successful
			{
				dt = dtBase;
				stepType = 0;
				meshId++;
			}
			if (stepType == 2) //if last step was a double mesh point
			{
				dt = 2 * meshPrecision;
				stepType = 0;
				meshId++;
			}
			if (nrOfEvents != 0 && stepType == 5) //dt should go back to normal
			{
				dt = dtBase;
				stepType = 0;
			}
			if (nrOfEvents != 0 && stepType == 3) //dt = dt / 2 and flush lastIndex
			{
				dt = 0.5 * dt;
				for (size_t i = 0; i < nrOfDelays; i++)
				{
					lastIndex[i] = 0;
				}
				if (dt < eventPrecision)
				{
					stepType = 5;
					nearEvent = false;
					dt = eventPrecision;

					//intervention in the integration
					for (size_t i = 0; i < nrOfUnroll; i++) //go through unroll
					{
						for (size_t j = 0; j < nrOfEvents; j++)
						{
							//check for positive direction of event
							Vec4d tmp = prevVals[nrOfEvents * i + j] * newVals[nrOfEvents * i + j];
							Vec4db mask = tmp < 0.0 && (prevVals[nrOfEvents * i + j] < newVals[nrOfEvents * i + j]);
							prevVals[nrOfEvents * i + j] = select(mask, NAN, prevVals[nrOfEvents * i + j]);
							int eventDir = 1;
							eventIntervention(j, eventDir, mask,t, x + nrOfVars * i, xDelay + nrOfDelays * i, p + nrOfParameters * i);

							//check for negativ direction of event
							mask = tmp < 0.0 && (prevVals[nrOfEvents * i + j] > newVals[nrOfEvents * i + j]);
							prevVals[nrOfEvents * i + j] = select(mask, NAN, prevVals[nrOfEvents * i + j]);
							eventDir = -1;
							eventIntervention(j, eventDir, mask,t, x + nrOfVars * i, xDelay + nrOfDelays * i, p + nrOfParameters * i);

						}
					}
				}
			}

			//detecting a mesh point
			if (meshId < meshLen && mesh[meshId] <= t + dt)
			{
				if(nrOfEvents != 0 && stepType == 5 ) //forced step
				{
					meshId++;
				}
				else
				{
					stepType = meshType[meshId];
					if (stepType == 1)
					{
						dt = mesh[meshId] - t;
					}
					if (stepType == 2)
					{
						dt = mesh[meshId] - t - meshPrecision;
					}
				}
			}//end of mesh detection

			//RK4 step
			RK4(f);

			//find event
			if (nrOfEvents != 0 && stepType != 5) //optimized out if there's no event
			{
				uint64_t sum = 0;

				//going through each variable to step for events
				for (size_t i = 0; i < nrOfUnroll; i++)
				{
					//push values, if step was not flushed
					if (stepType != 3)
					{
						for (size_t j = 0; j < nrOfEvents; j++)
						{
							prevVals[i * nrOfEvents + j] = newVals[i * nrOfEvents + j];
						}
					}
					eventLocation(newVals + nrOfEvents * i, t, xTmp + nrOfVars * i, xDelay + nrOfVars * i, p + nrOfParameters * i);
					//check for event
					for (size_t j = 0; j < nrOfEvents; j++)
					{
						Vec4d prod = prevVals[i * nrOfEvents + j] * newVals[i * nrOfEvents + j];
						Vec4db res = prod < 0; //1 when zero in the interval
						sum += horizontal_add(res);
					}
				}

				//sum > 0, step through event, unaccaptable step
				if (sum > 0)
				{
					stepType = 3;
					nearEvent = true;
				}
				//near an event, but we dont reach that yet, stepsize should stay
				if (sum == 0 && nearEvent == true)
				{
					stepType = 4;
				}
			} //end of find event

			//if step is acceptable
			if (nrOfEvents == 0 || stepType != 3)
			{

				t += dt;
				for (size_t i = 0; i < nrOfVars*nrOfUnroll; i++)
				{
					x[i] = xTmp[i];
				}

				//save end values
				memoryId++;
				tVals[memoryId] = t;
				for (size_t i = 0; i < nrOfUnroll; i++)
				{
					//save dense output points
					for (size_t j = 0; j < nrOfDenseVars; j++)
					{
						unsigned int varId = denseIdToVarId[j];
						x[i * nrOfVars + varId].store(xVals[memoryId][j] + i * vecSize);
					}
				}
			} //end of accept step
		} //end of while
	}//end of integrate

	void RK4(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		//k1
		calculateAllDelay(t);
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			f(&(kSum[nrOfVars * i]), t, &(x[nrOfVars * i]), &(xDelay[nrOfDelays * i]), &(p[nrOfParameters * i]));
		}
		//save dense output points
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				unsigned int varId = denseIdToVarId[j];
				kSum[i * nrOfVars + varId].store(xdVals[memoryId][j] + i * vecSize);
			}
		}


		//k2
		tTmp = t + 0.5 * dt;
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars * nrOfUnroll; i++)
		{
			xTmp[i] = x[i] + 0.5 * dt * kSum[i];
		}
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			f(&(kAct[nrOfVars * i]), tTmp, &(xTmp[nrOfVars * i]), &(xDelay[nrOfDelays * i]), &(p[nrOfParameters * i]));
		}

		//k3
		for (size_t i = 0; i < nrOfVars * nrOfUnroll; i++)
		{
			kSum[i] += 2 * kAct[i];
			xTmp[i] = x[i] + 0.5 * dt * kAct[i];
		}
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			f(&(kAct[nrOfVars * i]), tTmp, &(xTmp[nrOfVars * i]), &(xDelay[nrOfDelays * i]), &(p[nrOfParameters * i]));
		}


		//k4
		tTmp = t + dt;
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars * nrOfUnroll; i++)
		{
			kSum[i] += 2 * kAct[i];
			xTmp[i] = x[i] + dt * kAct[i];
		}
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			f(&(kAct[nrOfVars * i]), tTmp, &(xTmp[nrOfVars * i]), &(xDelay[nrOfDelays * i]), &(p[nrOfParameters * i]));
		}

		for (size_t i = 0; i < nrOfVars * nrOfUnroll; i++)
		{
			kSum[i] += kAct[i];
			xTmp[i] = x[i] + (1. / 6.) * dt * kSum[i];
		}
	}

	//dense output
	void calculatedenseOutput(double atT, int delayId)
	{
		//load new step if neccessary
		unsigned int id = lastIndex[delayId];
		unsigned int denseVar = delayIdToDenseId[delayId];
		if (atT >= tVals[id]) //new step id needed
		{
			while (atT >= tVals[id]) //find next good value
			{
				id++;
			}
			lastIndex[delayId] = id--;

			//load next step from memory
			tb[delayId] = tVals[id];
			deltat[delayId] = tVals[id + 1] - tb[delayId];
			pdeltat[delayId] = 1.0 / deltat[delayId];

			//load from memory
			for (size_t i = 0; i < nrOfUnroll; i++)
			{
				unsigned int linIndexOffset = nrOfDelays * i + delayId;
				xb[linIndexOffset].load(xVals[id][denseVar] + i * vecSize);
				xn[linIndexOffset].load(xVals[id + 1][denseVar] + i * vecSize);
				xdb[linIndexOffset].load(xdVals[id][denseVar] + i * vecSize);
				xdn[linIndexOffset].load(xdVals[id + 1][denseVar] + i * vecSize);
			}
		}

		//calculate dense output
		double theta = (atT - tb[delayId]) * pdeltat[delayId];
		double thetaM1 = theta - 1;
		for (size_t i = 0; i < nrOfUnroll; i++)
		{
			unsigned int linIndexOffset = nrOfDelays * i + delayId;
			xDelay[linIndexOffset] = -thetaM1 * xb[linIndexOffset] + theta * (xn[linIndexOffset] + thetaM1 * ((1 - 2 * theta) * (xn[linIndexOffset] - xb[linIndexOffset]) + deltat[delayId] * (thetaM1 * xdb[linIndexOffset] + theta * xdn[linIndexOffset])));
			}
	}
	void calculateAllDelay(double atT)
	{
		for (size_t i = 0; i < nrOfDelays; i++)
		{
			calculatedenseOutput(atT - t0[i], i);
		}
	}

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
		ofs.flush();
		ofs.close();
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
