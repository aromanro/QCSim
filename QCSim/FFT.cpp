#include "FFT.h"

namespace Fourier {
	const int FFT::init = fftw_init_threads();
	
	std::mutex FFTWPlan::planMutex;

	FFT::FFT(int numThreads)
	{		
		if (numThreads != 0) SetNumThreads(numThreads);
	}

	void FFT::SetNumThreads(int numThreads)
	{
		Clear();

		std::lock_guard lock(FFTWPlan::planMutex);

		fftw_plan_with_nthreads(numThreads);
	}

}