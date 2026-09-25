#pragma once

#include <exception>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace QC {

	// Runs func with the OpenMP parallel regions it starts from the calling thread limited to a single thread,
	// without affecting other threads. Eigen parallelizes the matrix products (for example the ones inside the SVDs)
	// with OpenMP, so this is what the MPS and MPO simulators use when their multithreading is disabled, which is useful
	// when several simulators run in parallel in different threads (otherwise each of them would start its own threads).
	//
	// The legacy MSVC OpenMP runtime (/openmp) keeps the number of threads setting for the whole process, so changing it
	// would affect the other threads as well, but it serializes the nested parallel regions, so func is run inside a
	// single thread parallel region.
	// The other runtimes (GCC, Clang, MSVC with /openmp:llvm) follow OpenMP >= 3.0: the number of threads is a per thread
	// setting, so it's set to one while func runs and restored afterwards (a nested parallel region would not be serialized there).
	//
	// NOTE: If Eigen::setNbThreads was called with a nonzero value, Eigen uses that instead of the OpenMP setting
	// (the legacy MSVC variant still works in that case, the other one doesn't).
	template<class Func> void RunSingleThreaded(Func&& func)
	{
#ifdef _OPENMP
#if defined(_MSC_VER) && !defined(__clang__) && !defined(_OPENMP_LLVM_RUNTIME)
		// an exception must not leave the parallel region, so it's passed out of it
		std::exception_ptr exception;

#pragma omp parallel num_threads(1)
		{
			try
			{
				func();
			}
			catch (...)
			{
				exception = std::current_exception();
			}
		}

		if (exception)
			std::rethrow_exception(exception);
#else
		struct ThreadsRestorer
		{
			const int previousThreads = omp_get_max_threads();
			~ThreadsRestorer() { omp_set_num_threads(previousThreads); }
		} restorer;

		omp_set_num_threads(1);

		func();
#endif
#else
		func();
#endif
	}

	// Runs func either normally or single threaded (see above)
	template<class Func> void RunMaybeSingleThreaded(bool multithreading, Func&& func)
	{
		if (multithreading)
			func();
		else
			RunSingleThreaded(func);
	}

}
