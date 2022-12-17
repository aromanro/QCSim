#include "FFT.h"

namespace Fourier {
	static int init = fftw_init_threads();
	
	std::mutex FFTWPlan::planMutex;

	FFT::FFT(int numThreads)
	{		
		if (numThreads != 0) SetNumThreads(numThreads);
	}


	FFT::~FFT()
	{
	}


	void FFT::SetNumThreads(int numThreads)
	{
		Clear();

		std::lock_guard<std::mutex> lock(FFTWPlan::planMutex);

		fftw_plan_with_nthreads(numThreads);
	}

}