#include "ParalellDDE.h"

template<unsigned int nrOfVars, unsigned int nrOfDelays, unsigned int nrOfDenseVars, unsigned int nrOfParameters, unsigned int nrOfEvents, unsigned int nrOfUnroll>
void ParalellDDE<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters, nrOfEvents, nrOfUnroll>::integrate(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
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
			x[i*nrOfVars + j] = x0[i * nrOfVars + j];
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
	meshId = 0;

	//integration
	int stepType;
	while(t < tEnd)
	{
		//assuming a simple step
		stepType = 0;
		dt = dtBase;

		//detecting a mesh point
		if (meshId < meshLen && mesh[meshId] <= t + dt)
		{
			stepType = meshType[meshId];
			std::cout << "step type: " << stepType << std::endl;
			if(stepType == 1) dt = mesh[meshId] - t;
			if (stepType == 2)
			{
				dt = mesh[meshId] - t - meshPrecision;
			}
			meshId++;
		}

		//RK4 step
		RK4step(f);


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
		//std::cout << std::setprecision(16) << "t = " << t << "\tx = " << x[0][0] << "\t" << x[1][0] << std::endl;

		if (stepType == 2)
		{
			std::cout << "Mini step" << std::endl;
			dt = 2 * meshPrecision;
			RK4step(f);


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
			std::cout << std::setprecision(16) << "t = " << t << "\tx = " << x[0][0] << "\t" << x[1][0] << std::endl;
		}

		
	}

}

template<unsigned int nrOfVars, unsigned int nrOfDelays, unsigned int nrOfDenseVars, unsigned int nrOfParameters, unsigned int nrOfEvents, unsigned int nrOfUnroll>
void ParalellDDE<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters, nrOfEvents, nrOfUnroll>::RK4step(void f(Vec4d* xd, double t, Vec4d* x, Vec4d* xDelay, Vec4d* p))
{
	//k1
	calculateAllDelay(t);
	for (size_t i = 0; i < nrOfUnroll; i++)
	{
		f(& (kSum[nrOfVars * i]), t, & (x[nrOfVars * i]), & (xDelay[nrOfDelays * i]), & (p[nrOfParameters * i]));
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


	//result of step
	for (size_t i = 0; i < nrOfVars * nrOfUnroll; i++)
	{
		kSum[i] += kAct[i];
		x[i] += (1. / 6.) * dt * kSum[i];
	}
	t += dt;

}

template<unsigned int nrOfVars, unsigned int nrOfDelays, unsigned int nrOfDenseVars, unsigned int nrOfParameters, unsigned int nrOfEvents, unsigned int nrOfUnroll>
inline void ParalellDDE<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters, nrOfEvents, nrOfUnroll>::calculatedenseOutput(double atT, int delayId)
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
			xb [linIndexOffset].load(xVals[id][denseVar] + i * vecSize);
			xn [linIndexOffset].load(xVals[id + 1][denseVar] + i * vecSize);
			xdb[linIndexOffset].load(xdVals[id][denseVar] + i * vecSize);
			xdn[linIndexOffset].load(xdVals[id + 1][denseVar] + i * vecSize);
			//std::cout << "Loaded tDelay = " << tDelay << " t_b= " << tb << " x_b= " << xb[i][0] << " x_n= " << xn[i][0] << " xd_b= " << xdb[i][0] << " xd_n= " << xdn[i][0] << std::endl;
		}
	}

	//calculate dense output
	double theta = (atT - tb[delayId]) * pdeltat[delayId];
	double thetaM1 = theta - 1;
	for (size_t i = 0; i < nrOfUnroll; i++)
	{
		unsigned int linIndexOffset = nrOfDelays * i + delayId;
		xDelay[linIndexOffset] = -thetaM1 * xb[linIndexOffset] + theta * (xn[linIndexOffset] + thetaM1 * ((1 - 2 * theta) * (xn[linIndexOffset] - xb[linIndexOffset]) + deltat[delayId] * (thetaM1 * xdb[linIndexOffset] + theta * xdn[linIndexOffset])));
		//std::cout << "th = " << theta << std::endl;
	}
}

template<unsigned int nrOfVars, unsigned int nrOfDelays, unsigned int nrOfDenseVars, unsigned int nrOfParameters, unsigned int nrOfEvents, unsigned int nrOfUnroll>
inline void ParalellDDE<nrOfVars, nrOfDelays, nrOfDenseVars, nrOfParameters, nrOfEvents, nrOfUnroll>::calculateAllDelay(double atT)
{
	for (size_t i = 0; i < nrOfDelays; i++)
	{
		calculatedenseOutput(atT - t0[i], i);
	}
}
