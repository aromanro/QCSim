#pragma once
#include <fftw3.h>

#include <complex>
#include <map>
#include <tuple>
#include <mutex>

// inspired from unsupported FFT from Eigen
// unfortunately the implementation there does not support 2D and 3D transforms and also multi threading
// I didn't like the way they index plans, either, so here it is, reimplemented

namespace Fourier {

	class FFTWPlan {
	public:
		FFTWPlan() = default;
		FFTWPlan(const FFTWPlan&) = delete;
		FFTWPlan& operator=(const FFTWPlan&) = delete;

		~FFTWPlan() 
		{ 
		  if (plan) 
		  { 
			std::lock_guard lock(planMutex);
			fftw_destroy_plan(plan); 
		  } 
		}

		// 1D

		inline void fwd(fftw_complex* src, fftw_complex* dst, unsigned int n)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_1d(n, src, dst, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void inv(fftw_complex* src, fftw_complex* dst, unsigned int n) {
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_1d(n, src, dst, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void fwd(double* src, fftw_complex* dst, unsigned int n)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_r2c_1d(n, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_r2c(plan, src, dst);
		}

		inline void inv(fftw_complex* src, double* dst, unsigned int n)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_c2r_1d(n, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_c2r(plan, src, dst);
		}

		// 2D

		inline void fwd(fftw_complex* src, fftw_complex* dst, unsigned int n0, unsigned int n1)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_2d(n0, n1, src, dst, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void inv(fftw_complex* src, fftw_complex* dst, unsigned int n0, unsigned int n1)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_2d(n0, n1, src, dst, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void fwd(double* src, fftw_complex* dst, unsigned int n0, unsigned int n1)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_r2c_2d(n0, n1, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_r2c(plan, src, dst);
		}

		inline void inv(fftw_complex* src, double* dst, unsigned int n0, unsigned int n1)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_c2r_2d(n0, n1, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_c2r(plan, src, dst);
		}

		// 3D

		inline void fwd(fftw_complex* src, fftw_complex* dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_3d(n0, n1, n2, src, dst, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void inv(fftw_complex* src, fftw_complex* dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			if (!plan) 
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_3d(n0, n1, n2, src, dst, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft(plan, src, dst);
		}

		inline void fwd(double* src, fftw_complex* dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			if (!plan) 
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_r2c_3d(n0, n1, n2, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_r2c(plan, src, dst);
		}

		inline void inv(fftw_complex* src, double* dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			if (!plan)
			{
				std::lock_guard lock(planMutex);
				plan = fftw_plan_dft_c2r_3d(n0, n1, n2, src, dst, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
			}
			fftw_execute_dft_c2r(plan, src, dst);
		}

		static std::mutex planMutex;

	private:
		fftw_plan plan = nullptr;
	};

	class FFT
	{
	public:
		explicit FFT(int numThreads = 0); // zero means let it alone for FFTW to decide

		// 1D

		inline void fwd(const std::complex<double>* src, std::complex<double> *dst, unsigned int n)
		{
			GetPlan(false, false, src, dst, n).fwd(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n);
		}

		inline void inv(const std::complex<double>* src, std::complex<double> *dst, unsigned int n)
		{
			GetPlan(true, false, src, dst, n).inv(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n);
		}

		inline void fwd(double* src, std::complex<double> *dst, unsigned int n)
		{
			GetPlan(false, true, src, dst, n).fwd(src, reinterpret_cast<fftw_complex*>(dst), n);
		}

		inline void inv(std::complex<double>* src, double *dst, unsigned int n)
		{
			GetPlan(true, true, src, dst, n).inv(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), dst, n);
		}

		// 2D

		inline void fwd(const std::complex<double>* src, std::complex<double> *dst, unsigned int n0, unsigned int n1)
		{
			GetPlan(false, false, src, dst, n0, n1).fwd(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n0, n1);
		}

		inline void inv(const std::complex<double>* src, std::complex<double> *dst, unsigned int n0, unsigned int n1)
		{
			GetPlan(true, false, src, dst, n0, n1).inv(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n0, n1);
		}


		// 3D

		inline void fwd(const std::complex<double>* src, std::complex<double> *dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			GetPlan(false, false, src, dst, n0, n1, n2).fwd(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n0, n1, n2);
		}

		inline void inv(const std::complex<double>* src, std::complex<double> *dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			GetPlan(true, false, src, dst, n0, n1, n2).inv(reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(src)), reinterpret_cast<fftw_complex*>(dst), n0, n1, n2);
		}

		void Clear()
		{
			Plans1D.clear();
			Plans2D.clear();
			Plans3D.clear();
			//fftw_cleanup_threads();
		}

		void SetNumThreads(int numThreads);

	private:
		// in place, aligned, inverse, different types, n
		std::map<std::tuple<bool, bool, bool, bool, unsigned int>, FFTWPlan> Plans1D;

		// in place, aligned, inverse, different types, n1, n2
		std::map<std::tuple<bool, bool, bool, bool, unsigned int, unsigned int>, FFTWPlan> Plans2D;

		// in place, aligned, inverse, different types, n1, n2, n3
		std::map<std::tuple<bool, bool, bool, bool, unsigned int, unsigned int, unsigned int>, FFTWPlan> Plans3D;

		inline bool InPlace(const void *src, const void *dst) { return src == dst; }
		inline bool Aligned(const void *src, const void *dst) { return ((reinterpret_cast<size_t>(src) & 0xF) | (reinterpret_cast<size_t>(dst) & 0xF)) == 0; }

		inline FFTWPlan& GetPlan(bool inverse, bool differentTypes, const void *src, void* dst, unsigned int n)
		{
			return Plans1D[std::tuple<bool, bool, bool, bool, unsigned int>(InPlace(src, dst), Aligned(src, dst), inverse, differentTypes, n)];
		}

		inline FFTWPlan& GetPlan(bool inverse, bool differentTypes, const void *src, void* dst, unsigned int n0, unsigned int n1)
		{
			return Plans2D[std::tuple<bool, bool, bool, bool, unsigned int, unsigned int>(InPlace(src, dst), Aligned(src, dst), inverse, differentTypes, n0, n1)];
		}

		inline FFTWPlan& GetPlan(bool inverse, bool differentTypes, const void *src, void* dst, unsigned int n0, unsigned int n1, unsigned int n2)
		{
			return Plans3D[std::tuple<bool, bool, bool, bool, unsigned int, unsigned int, unsigned int>(InPlace(src, dst), Aligned(src, dst), inverse, differentTypes, n0, n1, n2)];
		}

		static const int init;
	};


}