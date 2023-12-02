#pragma once

#include "QuantumAlgorithm.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <vector>
#include <stack>

namespace MachineLearning {

	// Implementation following the descriptin in Lesson 9 from:
	// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
	// https://arxiv.org/abs/2008.10647

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QMeansClustering2D :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		struct DataPoint {
			double x = 0.;
			double y = 0.;
		};

		struct ClusterDataPoint : DataPoint {
			int cluster = 0;
		};

		QMeansClustering2D(int addseed = 0)
		: BaseClass(3, addseed)
		{
			const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
			std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
			rng.seed(seed);
		}

		unsigned int Execute() override
		{
			ExecuteWithoutMeasurement();

			return BaseClass::Measure(0);
		}

		std::map<unsigned int, unsigned int> ExecuteWithMultipleMeasurements(unsigned int nrMeasurements = 1000)
		{
			ExecuteWithoutMeasurement();

			return BaseClass::RepeatedMeasure(0, 0, nrMeasurements);
		}

		void setData(const std::vector<DataPoint>& data)
		{
			m_data.resize(data.size());

			xmin = DBL_MAX;
			xmax = -DBL_MAX;
			ymin = DBL_MAX;
			ymax = -DBL_MAX;

			for (int i = 0; i < data.size(); ++i)
			{
				m_data[i].x = data[i].x;
				m_data[i].y = data[i].y;
				m_data[i].cluster = 0;

				xmin = std::min(xmin, data[i].x);
				xmax = std::max(xmax, data[i].x);
				ymin = std::min(ymin, data[i].y);
				ymax = std::max(ymax, data[i].y);
			}

			xspan = xmax - xmin;
			yspan = ymax - ymin;

			xspan = std::max(xspan, yspan);
			yspan = xspan;
		}

		const std::vector<ClusterDataPoint>& getData() const
		{
			return m_data;
		}

		const std::vector<ClusterDataPoint>& getCentroids() const
		{
			return centroids;
		}

		bool Cluster(unsigned int k, unsigned int nrMeasurements = 1000, unsigned int maxSteps = 100000)
		{
			getStartCentroids(k);

			return Update(nrMeasurements, maxSteps);
		}

		void SetDataPoint1(const ClusterDataPoint& dataPoint)
		{
			SetBlochStateParams(ugate1, dataPoint);
		}

		void SetDataPoint2(const ClusterDataPoint& dataPoint)
		{
			SetBlochStateParams(ugate2, dataPoint);
		}

	protected:
		void SetBlochStateParams(QC::Gates::UGate<MatrixClass>& ugate, const ClusterDataPoint& dataPoint)
		{
			const double x = (dataPoint.x - xmin) / xspan;
			const double y = (dataPoint.y - ymin) / yspan;

			const double theta = M_PI_2 * (x + y);
			const double phi = M_PI_2 * (x + 1. - y);

			ugate.SetParams(theta, phi);
		}

		void ExecuteWithoutMeasurement()
		{
			BaseClass::setToBasisState(0);

			BaseClass::ApplyGate(ugate1, 1);
			BaseClass::ApplyGate(ugate2, 2);

			// swap test
			BaseClass::ApplyGate(hadamard, 0);
			BaseClass::ApplyGate(cswap, 2, 1, 0);
			BaseClass::ApplyGate(hadamard, 0);
		}

		const std::vector<ClusterDataPoint>& getStartCentroids(unsigned int k)
		{
			static const double eps = 0.01 * xspan * xspan;
			centroids.resize(k);

			std::stack<unsigned int> s;

			unsigned int index = 0;
			bool retry = false;
			for (unsigned int c = 0; c < k; ++c)
			{
				std::uniform_int_distribution<> dist_nr(c, static_cast<unsigned int>(m_data.size()) - 1);
				
				do {
					index = dist_nr(rng);

					centroids[c].x = m_data[index].x;
					centroids[c].y = m_data[index].y;
					centroids[c].cluster = c;

					// don't allow two centroids to be too close to each other
					retry = false;
					for (unsigned int oc = 0; oc < c; ++oc)
					{
						const double dx = centroids[c].x - centroids[oc].x;
						const double dy = centroids[c].y - centroids[oc].y;

						if (dx * dx + dy * dy < eps)
						{
							retry = true;
							break;
						}
					}
				} while (retry);

				s.push(index);

				// move the selected data point to the beginning, to avoid selecting it again
				std::swap(m_data[c], m_data[index]);
			}

			// restore the order of data points, just in case someone expects them to remain in the same order
			while (!s.empty())
			{
				std::swap(m_data[s.top()], m_data[s.size() - 1]);
				s.pop();
			}

			return centroids;
		}

		void UpdateCentroids()
		{
			std::vector<int> clusterSizes(centroids.size(), 0);

			for (int c = 0; c < centroids.size(); ++c)
				centroids[c].x = centroids[c].y = 0;

			for (const auto& p : m_data)
			{
				centroids[p.cluster].x += p.x;
				centroids[p.cluster].y += p.y;
				++clusterSizes[p.cluster];
			}

			for (int c = 0; c < centroids.size(); ++c)
			{
				if (clusterSizes[c] == 0) continue;
				centroids[c].x /= clusterSizes[c];
				centroids[c].y /= clusterSizes[c];
			}
		}

		bool UpdateClusters(unsigned int nrMeasurements)
		{
			bool terminate = true;

			for (int i = 0; i < m_data.size(); ++i)
			{
				unsigned int overlap = 0;
				int cluster = 0;

				SetDataPoint1(m_data[i]);

				for (int c = 0; c < centroids.size(); ++c)
				{
					SetDataPoint2(centroids[c]);
				
					auto res = ExecuteWithMultipleMeasurements(nrMeasurements);
					assert(res.size() == 1 || res.size() == 2);

					if (res[0] > overlap)
					{
						overlap = res[0];
						cluster = c;
					}
				}

				if (m_data[i].cluster != cluster)
				{
					terminate = false;
					m_data[i].cluster = cluster;
				}
			}

			return terminate;
		}

		bool Update(unsigned int nrMeasurements, unsigned int maxSteps)
		{
			bool terminate = false;

			for (unsigned int i = 0; i < maxSteps && !terminate; ++i)
				terminate = UpdateStep(nrMeasurements);

			return terminate;
		}

		bool UpdateStep(unsigned int nrMeasurements)
		{
			const bool terminate = UpdateClusters(nrMeasurements);
			if (!terminate) UpdateCentroids();

			return terminate;
		}

		std::vector<ClusterDataPoint> m_data;
		std::vector<ClusterDataPoint> centroids;

		double xmin = 0.;
		double xmax = 0.;
		double ymin = 0.;
		double ymax = 0.;
		double xspan = 1.;
		double yspan = 1.;

		std::mt19937_64 rng;

		QC::Gates::UGate<MatrixClass> ugate1;
		QC::Gates::UGate<MatrixClass> ugate2;
		QC::Gates::HadamardGate<MatrixClass> hadamard;
		QC::Gates::FredkinGate<MatrixClass> cswap;
	};

}


