#ifndef PARALLELDDE
#define PARALLELDDE

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "C:/Users/nnagy/Documents/Egyetem/HDS/VCL/version2-2.01.02/vectorclass.h"
#include "C:/Users/nnagy/Documents/Egyetem/HDS/VCL/version2-2.01.02/vectormath_exp.h"
#include "ParalellDDE_Common.cpp"
#include "DP5_Constants.h"
#define vecSize 4

template<unsigned int nrOfVars, unsigned int nrOfDelays = 0, unsigned int nrOfDenseVars = 0, unsigned int nrOfParameters = 0, unsigned int nrOfEvents = 0, unsigned int nrOfUnroll = 1>
class ParalellDDE_DP5
{
private:
	//sizes
	unsigned int nrOfSteps;
	unsigned int memoryId;

	//integration memory to every variable with dense output
	double*** __restrict denseK1, *** __restrict denseK2, *** __restrict denseK3, *** __restrict denseK4, *** __restrict denseK5;
	double***__restrict denseK6, *** __restrict denseX, ** __restrict denseT;
	/*
	Integration memory layout
	length of outer dimension:	denseOutputMemorySize
	length of middle dimension: nrOfDenseVars
	length of inner dimension:	vecSize
	*/

	//mesh points -> synchronized
	unsigned int meshLen; //number of mesh points
	unsigned int * meshId; //index of next mesh point
	double* mesh; //synchronized mesh points for variables -> avoid divergence
	Vec4d nextMesh;
	double meshPrecision;
	int* meshType; //1 = C1 disc or lower (two points saved)	2 = C2 disc. or higher (one point saved, stepped on that point)
	/*
	Every equation tracks it's mesh points
	*/


	//delays -> synchronized
	double * t0; 
	unsigned int delayIdToDenseId[nrOfDelays]; //	lookup table:	delay index		-> dense variable index
	unsigned int varIdToDenseId[nrOfVars]; //		lookup table:	variable index	-> dense variable index
	unsigned int denseIdToVarId[nrOfDenseVars]; //	lookup table:	dense variable	-> variable index
	unsigned int delayIdToVarId[nrOfDenseVars]; //	lookup table:	delay index	-> variable index

	//dynamic integration variables
	Vec4d *__restrict x;
	Vec4d t, dt; //synchronized for all threads

	//initial condition
	Vec4d* x0;
	Vec4d* endvals;
	Vec4d stepCounter;

	//temporary integration variables
	Vec4d *__restrict k1, * __restrict k2, * __restrict k3, * __restrict k4, * __restrict k5, * __restrict k6, * __restrict k7;
	Vec4d *__restrict xTmp, * __restrict xDelay, * __restrict x4, * __restrict x5, tTmp;
	/*
	integration variable layout
	indexing: unroll*nrOfvars + var
	*/

	//dense output memory -> to every delay
	Vec4d * __restrict lk1, * __restrict lk2, * __restrict lk3, * __restrict lk4, * __restrict lk5, * __restrict lk6, * __restrict theta, * __restrict lX, * __restrict lDt;
	unsigned int * __restrict lastIndex, * __restrict newLastIndex;

	//user given integration range
	double tStart, tEnd, dtBase;

	//parameters
	Vec4d * __restrict p;

	//initial discontinouities
	double* discC0Init, *discC1Init, *discInit;
	unsigned int nrOfC0, nrOfC1, nrOfDisc;

	//initial functions
	typedef double(*funcSlot)(double, int);
	funcSlot* initialFunction;

	//event stuff
	Vec4d* prevVals;
	Vec4d* newVals;
	double eventPrecision;
	void (*eventLocation)(Vec4d* lst, Vec4d t, Vec4d* x, Vec4d* xDelay, Vec4d* p);
	void (*eventIntervention)(int eventId, int eventDir, Vec4db mask, Vec4d t, Vec4d* x, Vec4d* xDelay, Vec4d* p);

	//adaptive solver
	double absTol, relTol;
	double dtMin, dtMax;
	double boundIncrease, boundDecrease;
	double safetyFactor;

public:
	//constructor
	ParalellDDE_DP5() : meshPrecision(1e-14), meshLen(0), nrOfC0(0), nrOfC1(0), nrOfDisc(0), nrOfSteps(100), dtBase(0.2), tStart(0.0), tEnd(10.0), eventPrecision(1e-10), absTol(1e-6), relTol(1e-6), dtMin(1e-10), dtMax(1e6), boundDecrease(0.1), boundIncrease(10.0), safetyFactor(0.8)
	{
		//initial function
		initialFunction = new funcSlot[nrOfDenseVars];

		//index prediction
		lastIndex = new unsigned int[nrOfDelays * vecSize];
		newLastIndex = new unsigned int[nrOfDelays * vecSize];
		meshId = new unsigned int[vecSize];
		for (size_t i = 0; i < nrOfDelays * vecSize; i++)
		{
			lastIndex[i] = 0;
		}

		//events
		prevVals = new Vec4d[nrOfEvents];
		newVals = new Vec4d[nrOfEvents];

		//integration
		p = new Vec4d[nrOfParameters];
		x = new Vec4d[nrOfVars];
		endvals = new Vec4d[nrOfVars];
		x0 = new Vec4d[nrOfVars];
		x4 = new Vec4d[nrOfVars];
		x5 = new Vec4d[nrOfVars];
		k1 = new Vec4d[nrOfVars];
		k2 = new Vec4d[nrOfVars];
		k3 = new Vec4d[nrOfVars];
		k4 = new Vec4d[nrOfVars];
		k5 = new Vec4d[nrOfVars];
		k6 = new Vec4d[nrOfVars];
		k7 = new Vec4d[nrOfVars];
		xTmp = new Vec4d[nrOfVars];
		xDelay = new Vec4d[nrOfDelays];
		
		//dense output
		lk1 =   new Vec4d[nrOfDelays];
		lk2 =   new Vec4d[nrOfDelays];
		lk3 =   new Vec4d[nrOfDelays];
		lk4 =   new Vec4d[nrOfDelays];
		lk5 =   new Vec4d[nrOfDelays];
		lk6 =   new Vec4d[nrOfDelays];
		lX =    new Vec4d[nrOfDelays];
		lDt =   new Vec4d[nrOfDelays];
		theta = new Vec4d[nrOfDelays];

		//delays
		t0 = new double[nrOfDelays];
	};

	//destruktor
	~ParalellDDE_DP5()
	{
		/*delete p, x, x0, kAct, kSum, xTmp, xDelay;
		delete xb, xn, xdb, xdn, tb, deltat, pdeltat, t0;*/
	}; //does nothing

	//set functions
	void setSolverTol(double absTol, double relTol)
	{
		this->absTol = absTol;
		this->relTol = relTol;
	}
	void setGlobalTol(double tol)
	{
		this->absTol = tol;
		this->relTol = tol;
		this->meshPrecision = tol;
		this->eventPrecision = tol;
	}
	void setStepSizeBounds(double dtMin, double dtMax)
	{
		this->dtMin = dtMin;
		this->dtMax = dtMax;
	}
	void setStepSizeChangeBounds(double decrease, double increase)
	{
		this->boundIncrease = increase;
		this->boundDecrease = decrease;
	}
	void setSafetyFactor(double safetyFactor)
	{
		this->safetyFactor = safetyFactor;
	}
	void setIntegrationRange(double tStart, double tEnd)
	{
		this->tStart = tStart;
		this->tEnd = tEnd;
	}
	void setStepSize(double dt)
	{
		this->dtBase = dt;
	}
	void setNrOfSteps(unsigned int nrOfSteps)
	{
		this->nrOfSteps = nrOfSteps;
		allocateMemory();
	}
	void setNrOfDisc(unsigned int nrOfDisc)
	{
		this->nrOfDisc = nrOfDisc;
	}
	void allocateMemory()
	{
		//allocate arrays
		this->denseT = new double* [nrOfSteps];
		this->denseX = new double** [nrOfSteps];
		this->denseK1 = new double** [nrOfSteps];
		this->denseK2 = new double** [nrOfSteps];
		this->denseK3 = new double** [nrOfSteps];
		this->denseK4 = new double** [nrOfSteps];
		this->denseK5 = new double** [nrOfSteps];
		this->denseK6 = new double ** [nrOfSteps];

		for (size_t i = 0; i < nrOfSteps; i++)
		{
			this->denseT [i] = new double[vecSize];
			for (size_t j = 0; j < vecSize; j++)
			{
				this->denseT[i][j] = NAN;
			}
			this->denseX [i] = new double* [nrOfDenseVars];
			this->denseK1[i] = new double* [nrOfDenseVars];
			this->denseK2[i] = new double* [nrOfDenseVars];
			this->denseK3[i] = new double* [nrOfDenseVars];
			this->denseK4[i] = new double* [nrOfDenseVars];
			this->denseK5[i] = new double* [nrOfDenseVars];
			this->denseK6[i] = new double* [nrOfDenseVars];
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				this->denseX [i][j] = new double[vecSize];
				this->denseK1[i][j] = new double[vecSize];
				this->denseK2[i][j] = new double[vecSize];
				this->denseK3[i][j] = new double[vecSize];
				this->denseK4[i][j] = new double[vecSize];
				this->denseK5[i][j] = new double[vecSize];
				this->denseK6[i][j] = new double[vecSize];
				for (size_t k = 0; k < vecSize; k++)
				{
					this->denseX[i][j] [k]=	NAN;
					this->denseK1[i][j][k] =NAN;
					this->denseK2[i][j][k] =NAN;
					this->denseK3[i][j][k] =NAN;
					this->denseK4[i][j][k] =NAN;
					this->denseK5[i][j][k] =NAN;
					this->denseK6[i][j][k] =NAN;
				}
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
		this->delayIdToDenseId[delayId] = denseVarId; //delay --> dense var
		this->delayIdToVarId[delayId] = this->denseIdToVarId[denseVarId]; //delay --> var
	}
	void addInitialFunction(int varId, int denseVarId, double f(double, int))
	{
		this->initialFunction[denseVarId] = f;
		this->varIdToDenseId[varId] = denseVarId; //var --> dense
		this->denseIdToVarId[denseVarId] = varId; //dense --> var

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
	unsigned int getDenseOutputMemorySize() //returns whole memory size in bytes
	{
		return 8 * nrOfSteps * nrOfDenseVars  * vecSize * sizeof(double);
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
		double* endVals = new double[vecSize];

		for (size_t j = 0; j < vecSize; j++)
		{
			endVals[j] = this->endvals[varId][j];
		}
		return endVals;
	}
	double getT()
	{
		return t;
	}
	int* getStepCount()
	{
		int* stepCount = new int[vecSize];
		for (size_t i = 0; i < vecSize; i++)
		{
			stepCount[i] = int(this->stepCounter[i]);
		}
		return stepCount;
	}

	//calculate functions - alternative to set, it does the job itself
	void calculateIntegrationMesh()
	{
		double* discC0Full = new double[nrOfC0 + 1];
		discC0Full[0] = 0.0;
		for (size_t i = 1; i < nrOfC0 + 1; i++)
		{
			discC0Full[i] = discC0Init[i - 1];
		}
		calculateFullIntegrationMesh(&mesh, &meshType, &meshLen, discC0Full, nrOfC0 + 1, discC1Init, nrOfC1, t0, nrOfDelays, tStart, meshPrecision);
		this->addMeshPoint(tEnd, 1);
	}

	//integration functions
	void integrate(void f(Vec4d* xd, Vec4d t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		// ---------------------------- Set everything to initial values -----------------------------------------
		//initialize variables
		t = tStart; //tStart -> t
		memoryId = 0;
		dt = dtBase;
		for (size_t i = 0; i < nrOfVars; i++)
		{
			x[i] = x0[i];
			endvals[i] = 0;
		}
		if (meshLen > 0)
		{
			nextMesh = mesh[0];
		}

		//save dense output points
		t.store(denseT[memoryId]);
		for (size_t i = 0; i < nrOfDenseVars; i++)
		{
			unsigned int varId = denseIdToVarId[i];
			x[i].store(denseX[memoryId][i]);
		}
		for (size_t i = 0; i < nrOfDelays * vecSize; i++) //delay index prediction
		{
			lastIndex[i] = 0;
			newLastIndex[i] = 0;
		}
		for (size_t i = 0; i < nrOfEvents; i++) //event handling
		{
			prevVals[i] = NAN;
			newVals[i] = NAN;
		}
		for (size_t i = 0; i < vecSize; i++) //mesh
		{
			meshId[i] = 0;
		}
		
		//integration
		//0 normal step
		//1 simple mesh point
		//2 double mesh point
		//3 stepsize decrease because event, decline step
		//4 acceptable step near event
		//5 event is found, we are after that
		//6 decline step
		//7 dt should reset
		Vec4d stepType = 0;
		Vec4d nearEvent = 0;
		Vec4d newDt = dtBase;
		Vec4d saved = 0;
		for(memoryId = 1; memoryId < nrOfSteps; memoryId++)
		{
			for (size_t i = 0; i < vecSize; i++)
			{
				if (stepType[i] == 1 || stepType[i] == 2)
				{
					meshId[i]++;
				}
			}

			dt = select(stepType == 0 | stepType == 6, newDt, dt); //when dt stays
			dt = select(stepType == 1 | stepType == 5 | stepType == 7, dtBase, dt); //reset dt
			dt = select(stepType == 2, 2 * meshPrecision, dt);
			dt = select(stepType == 3, 0.5 * dt, dt);

			stepType = select(stepType == 1 | stepType == 5 | stepType == 6, 0, stepType);
			stepType = select(stepType == 2, 7, stepType);
			
			
			//event intervention
			Vec4db eventI = stepType == 3 & dt < eventPrecision;

			if (horizontal_add(eventI)) //if there is an event
			{
				stepType = select(eventI, 5, stepType);

				//intervention in the integration
				for (size_t j = 0; j < nrOfEvents; j++)
				{
					//check for positive direction of event
					Vec4d tmp = prevVals[j] * newVals[j];
					Vec4db mask = tmp < 0.0 && (prevVals[j] < newVals[j]);
					prevVals[j] = select(mask, NAN, prevVals[j]);
					int eventDir = 1;
					eventIntervention(j, eventDir, mask, t, x, xDelay, p);

					//check for negativ direction of event
					mask = tmp < 0.0 && (prevVals[j] > newVals[j]);
					prevVals[j] = select(mask, NAN, prevVals[j]);
					eventDir = -1;
					eventIntervention(j, eventDir, mask, t, x, xDelay, p);

				}
			}

			//detecting a mesh point this part is not parallel :((((((((
			for (size_t i = 0; i < vecSize; i++)
			{
				unsigned int id = meshId[i];
				if (id < meshLen && mesh[id] <= t[i] + dt[i]) 
				{
					if (stepType[i] == 5) //forced step
					{
						meshId[i]++;
					}
					else
					{

						int type = meshType[id];
						if (type == 1)
						{
							stepType.insert(i, 1);
							dt.insert(i, mesh[id] - t[i]);
						}
						if (type == 2)
						{
							stepType.insert(i, 2);
							dt.insert(i, mesh[id] - t[i] - meshPrecision);
						}
					}
				}
			}

			//RK4 step
			DP5(f);

			//find event
			if (nrOfEvents != 0) //optimized out if there's no event
			{
				uint64_t sum = 0;

				//push values, if step was not flushed
				for (size_t j = 0; j < nrOfEvents; j++)
				{
					prevVals[j] = select(stepType != 3 & stepType != 5, newVals[j], prevVals[j]);
				}
				eventLocation(newVals, tTmp, x5, xDelay, p);

				//check for event
				Vec4db whereEvent = false;
				for (size_t j = 0; j < nrOfEvents; j++)
				{
					Vec4d prod = prevVals[j] * newVals[j];
					Vec4db res = prod < 0; //1 when zero in the interval
					whereEvent = whereEvent | res;
				}

				//if a variable reached an event
				nearEvent = select(whereEvent & stepType != 5., 1, nearEvent);
				stepType = select(whereEvent & stepType != 5, 3, stepType);
				stepType = select(nearEvent == 1 & !whereEvent & stepType != 5, 4, stepType);
			} //end of find event

			//calculate new dt
			Vec4d minRelativeErrorReciprok = boundIncrease;
			for (size_t i = 0; i < nrOfVars; i++)
			{
				Vec4d errori = select(x4[i] < x5[i],x5[i] - x4[i], x4[i]-x5[i]);
				Vec4d tolerancei = absTol + select(x5[i]<0,-x5[i],x5[i]) * relTol;
				minRelativeErrorReciprok = min(minRelativeErrorReciprok, tolerancei / errori);
			}
			Vec4d changeFactor = 0.8*pow<double>(minRelativeErrorReciprok, 0.2);
			//check change factor bounds
			changeFactor = select(changeFactor < boundDecrease, boundDecrease, changeFactor);
			changeFactor = select(changeFactor > boundIncrease, boundIncrease, changeFactor);
			//new dt
			newDt = changeFactor * dt;
			//check dt bounds
			newDt = select(newDt < dtMin, dtMin, newDt);
			newDt = select(newDt > dtMax, dtMax, newDt);

			//check if step should be rejected
			stepType = select(minRelativeErrorReciprok < 1.0, 6, stepType);
			//check if dt is already small
			stepType = select(dt <= dtMin, 0, stepType);

			//accept or decline step
			Vec4db acceptStep = (stepType != 3 & stepType != 6) | stepType == 2;
			t = select(acceptStep, tTmp, t);
			for (size_t i = 0; i < nrOfVars; i++)
			{
				x[i] = select(acceptStep, x5[i], x[i]);
			}

			//if step accepted load new 
			for (size_t i = 0; i < nrOfDelays; i++)
			{
				for (size_t j = 0; j < vecSize; j++)
				{
					unsigned int linearIndex = i * vecSize + j;
					if (acceptStep[j] == 1 && newLastIndex[linearIndex] > 1)
					{
						lastIndex[linearIndex] = newLastIndex[linearIndex] - 1;
					}
				}
			}
			

			//save dense values
			t.store(denseT[memoryId]);
			for (size_t j = 0; j < nrOfDenseVars; j++)
			{
				unsigned int varId = denseIdToVarId[j];
				x [varId].store(denseX [memoryId][j]);
				k1[varId].store(denseK1[memoryId-1][j]);
				k2[varId].store(denseK2[memoryId-1][j]);
				k3[varId].store(denseK3[memoryId-1][j]);
				k4[varId].store(denseK4[memoryId-1][j]);
				k5[varId].store(denseK5[memoryId-1][j]);
				k6[varId].store(denseK6[memoryId-1][j]);
			}

			/*std::cout << "0:  t=" << t[0] << "\tx=" << x[0][0] << "\tdt=" << dt[0] << "\t"; 
			std::cout << "1:  t=" << t[1] << "\tx=" << x[0][1] << "\tdt=" << dt[1] << "\t"; 
			std::cout << "2:  t=" << t[2] << "\tx=" << x[0][2] << "\tdt=" << dt[2] << "\t"; 
			std::cout << "3:  t=" << t[3] << "\tx=" << x[0][3] << "\tdt=" << dt[3] << "\n";*/
			//std::cout << "t=" << t[1] << "\tx=" << x[0][1] << "\tdt=" << dt[1] << "\tk1="<<k1[0][1]<< "\tstep=" << stepType[0] << "\tError: "<<minRelativeErrorReciprok[0]<< "\n";
			/*std::cout << "\tk1= " << k1[0][0];
			std::cout << "\tk2= " << k2[0][0];
			std::cout << "\tk3= " << k3[0][0];
			std::cout << "\tk4= " << k4[0][0];
			std::cout << "\tk5= " << k5[0][0];
			std::cout << "\tk6= " << k6[0][0];
			std::cout << "\tid= " << memoryId;
			std::cout << std::endl;*/

			//check end
			if (horizontal_add(t >= tEnd))
			{
				for (size_t i = 0; i < nrOfVars; i++)
				{
					endvals[i] = select(saved == 1, endvals[i], x[i]);
					//std::cout << endvals[i][0] << std::endl;
				}
				stepCounter = select(saved == 1, stepCounter, double(memoryId));
				saved = select(t >= tEnd, 1, 0);
				/*std::cout << saved[0] << "\t" << saved[1] << "\t" << saved[2] << "\t" << saved[3] << "\t" << "\n";
				std::cout << stepCounter[0] << "\t" << stepCounter[1] << "\t" << stepCounter[2] << "\t" << stepCounter[3] << "\t" << "\n";
				std::cout << endvals[0][0] << "\t" << endvals[0][1] << "\t" << endvals[0][2] << "\t" << endvals[0][3] << "\t" << "\n";*/
				if (horizontal_and(t >= tEnd))
				{
					//std::cout << t[0] << "\t" << t[1] << "\t" << t[2] << "\t" << t[3] << "\t";
					break;
				}
			}

		} //end of while
	}//end of integrate
	void DP5(void f(Vec4d* xd, Vec4d t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
	{
		//k1
		calculateAllDelay(t);
		f(k1, t, x, xDelay, p);


		//k2
		tTmp = t + c1 * dt;	//calculate new t
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (a11 * k1[i]);
		}
		f(k2, tTmp, xTmp, xDelay, p);
		

		//k3
		tTmp = t + c2 * dt;	//calculate new t
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (a21 * k1[i] + a22 * k2[i]);
		}
		f(k3, tTmp, xTmp, xDelay, p);

		//k4
		tTmp = t + c3 * dt;	//calculate new t
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (a31 * k1[i] + a32 * k2[i] + a33 * k3[i]);
		}
		f(k4, tTmp, xTmp, xDelay, p);

		//k5
		tTmp = t + c4 * dt;	//calculate new t
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i] + a44 * k4[i]);
		}
		f(k5, tTmp, xTmp, xDelay, p);

		//k6
		tTmp = t + c5 * dt;	//calculate new t
		calculateAllDelay(tTmp);
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i] + a55 * k5[i]);
		}
		f(k6, tTmp, xTmp, xDelay, p);

		//k7 === x5
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			xTmp[i] = x[i] + dt * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i] + b55 * k5[i] + b56 * k6[i]);
			x5[i] = xTmp[i];
		}
		f(k7, tTmp, xTmp, xDelay, p);

		//x4
		for (size_t i = 0; i < nrOfVars; i++) //calculate new x values
		{
			x4[i] = x[i] + dt * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i] + b44 * k4[i] + b45 * k5[i] + b46 * k6[i] + b47 * k7[i]);
		}
	}
	
	//dense output
	void denseOutput(Vec4db mask, int delayId) //from values loaded into memory
	{
		Vec4d b1 = theta[delayId] * (1 + theta[delayId] * ((-1337.0 / 480.0) + theta[delayId] * ((1039.0 / 360.0) - theta[delayId] * (1163.0 / 1152.0))));
		Vec4d b2 = 0;
		Vec4d b3 = theta[delayId] * theta[delayId] * ((4216.0 / 1113.0) + theta[delayId] * ((-18728.0 / 3339.0) + theta[delayId] * (7580.0 / 3339.0)));
		Vec4d b4 = theta[delayId] * theta[delayId] * ((-27.0 / 16.0) + theta[delayId] * (( 9.0 / 2.0) - theta[delayId] * (415.0 / 192.0)));
		Vec4d b5 = theta[delayId] * theta[delayId] * ((-2187.0 / 8480.0) + theta[delayId] * ((2673.0 / 2120.0) - theta[delayId] * (8991.0 / 6784.0)));
		Vec4d b6 = theta[delayId] * theta[delayId] * ((33.0 / 35.0) + theta[delayId] * ((-319.0 / 105.0) + theta[delayId] * (187.0 / 84.0)));
		Vec4d tmp = lX[delayId] + lDt[delayId] * (b1 * lk1[delayId] + b2 * lk2[delayId] + b3 * lk3[delayId] + b4 * lk4[delayId] + b5 * lk5[delayId] + b6 * lk6[delayId]);
		xDelay[delayId] = select(mask,tmp , xDelay[delayId]);
	}
	void loadDenseOutput(double atT, int delayId, int vecId)
	{
		unsigned int linearIndex = delayId * vecSize + vecId;
		unsigned int id = lastIndex[linearIndex];
		unsigned int denseVar = delayIdToDenseId[delayId];

		if (atT >= denseT[id][vecId])  //new values should be loaded
		{
			while (atT >= denseT[id][vecId])//find next good value
			{
				id++;
			}
			newLastIndex[linearIndex] = id--;

			//load next step from memory
			lk1[delayId].insert(vecId, denseK1[id][denseVar][vecId]);
			lk2[delayId].insert(vecId, denseK2[id][denseVar][vecId]);
			lk3[delayId].insert(vecId, denseK3[id][denseVar][vecId]);
			lk4[delayId].insert(vecId, denseK4[id][denseVar][vecId]);
			lk5[delayId].insert(vecId, denseK5[id][denseVar][vecId]);
			lk6[delayId].insert(vecId, denseK6[id][denseVar][vecId]);
			lX[delayId].insert(vecId, denseX[id][denseVar][vecId]);

			//calculate theta
			double tBefore = denseT[id][vecId];
			double dtStep = denseT[id + 1][vecId] - tBefore;
			theta[delayId].insert(vecId, (atT - tBefore) / dtStep);
			lDt[delayId].insert(vecId, dtStep);
		}

	}
	void calculateAllDelay(Vec4d atT)
	{
		for (size_t i = 0; i < nrOfDelays; i++)
		{
			//going through variables, this part is not vectorized :(
			unsigned int denseId = delayIdToDenseId[i];
			unsigned int varId = delayIdToVarId[i];
			Vec4db mask;
			for (size_t j = 0; j < vecSize; j++)
			{
				Vec4d getT = atT - t0[i];
				loadDenseOutput(getT[j], i, j);
				//std::cout << "atT=" << getT[j] << " Theta=" << theta[0][0] << " x=" << lX[0][0] << std::endl;
				if (getT[j] <= tStart) //initial function
				{
					mask.insert(j, false);
					xDelay[i].insert(j, initialFunction[denseId](getT[j],0));
				}
				else //dense output
				{
					mask.insert(j, true);
				}
			}
			//std::cout << mask[0] << mask[1] << mask[2] << mask[3] << "\n";
			denseOutput(mask,i);
		}
		//std::cout << "t=" << atT[0] << " xd= "<< xDelay[0][1] <<  std::endl;
	}

	//save functions
	void save(std::string filename, unsigned int varId = 0, unsigned int vecId = 0, int precision = 16)
	{
		std::ofstream ofs(filename);
		if (!ofs.is_open())
		{
			std::cout << "Cannot open file: " << filename << std::endl;
			return;
		}
		ofs << std::setprecision(precision);
		for (size_t i = 0; i < memoryId; i++)
		{
			ofs << denseT[i][vecId] << "\t";
			ofs << denseX[i][varId][vecId] << std::endl;
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
