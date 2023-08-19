#include "Tests.h"

#include "QMeansClustering2D.h"

#include <algorithm>

int clusterForPoint(const MachineLearning::QMeansClustering2D<>::DataPoint& p, const std::vector<MachineLearning::QMeansClustering2D<>::DataPoint>& centroids)
{
	int res = -1;

	double dist = std::numeric_limits<double>::max();

	for (int i = 0; i < centroids.size(); ++i)
	{
		const double dx = p.x - centroids[i].x;
		const double dy = p.y - centroids[i].y;
		const double d2 = dx * dx + dy * dy;

		if (d2 < dist)
		{
			dist = d2;
			res = i;
		}
	}

	return res;
}

bool QMeansClustering2DTests()
{
	std::cout << "\nTesting Q-Means clustering..." << std::endl;

	const unsigned int nrMeasurements = 100000;

	unsigned int k = 3;
	unsigned int pointsPerCluster = 15;

	std::vector<MachineLearning::QMeansClustering2D<>::DataPoint> origCentroids(k);
	origCentroids[0].x = 0;
	origCentroids[0].y = 0.;

	origCentroids[1].x = 3;
	origCentroids[1].y = 3;

	origCentroids[2].x = 1.5;
	origCentroids[2].y = 2.;

	std::vector<MachineLearning::QMeansClustering2D<>::DataPoint> data;
	data.reserve(k * pointsPerCluster);

	std::uniform_real_distribution<> dist(-1.5, 1.5);

	for (unsigned int c = 0; c < k; ++c)
	{
		for (unsigned int j = 0; j < pointsPerCluster; ++j)
		{
			MachineLearning::QMeansClustering2D<>::DataPoint p;
			p.x = origCentroids[c].x + dist(gen);
			p.y = origCentroids[c].y + dist(gen);
			
			data.push_back(p);
		}
	}

	std::shuffle(data.begin(), data.end(), gen);
	
	std::cout << "Running clustering..." << std::endl;

	MachineLearning::QMeansClustering2D<> qMeansClustering;
	qMeansClustering.setData(data);
	if (!qMeansClustering.Cluster(k, nrMeasurements)) std::cout << "Did not terminate" << std::endl;

	std::cout << "Checking results..." << std::endl;

	std::vector<int> clusterCounts(k, 0);
	const auto& dataPoints = qMeansClustering.getData();
	for (int i = 0; i < dataPoints.size(); ++i)
		++clusterCounts[dataPoints[i].cluster];

	const auto& centroids = qMeansClustering.getCentroids();

	for (int c = 0; c < centroids.size(); ++c)
		std::cout << "Cluster " << c << " centroid: x: " << centroids[c].x << ", y: " << centroids[c].y << " Points count: " << clusterCounts[c] << std::endl;

	// not surprisingly it fails sometimes (especially for k > 2) for the euclidian metric test (since it's not the same metric as the one used for clustering)
	// but also sometimes (not so often as for the euclidian metric check, though) for the 'overlap' one, due of the probabilistic nature
	// increasing the number of measurements helps, of course, but obviously it slows down the execution
	
	for (int i = 0; i < dataPoints.size(); ++i)
	{
		const double dx = dataPoints[i].x - centroids[dataPoints[i].cluster].x;
		const double dy = dataPoints[i].y - centroids[dataPoints[i].cluster].y;
		const double d2 = dx * dx + dy * dy;

		qMeansClustering.SetDataPoint1(dataPoints[i]);
		qMeansClustering.SetDataPoint2(centroids[dataPoints[i].cluster]);
		auto res = qMeansClustering.ExecuteWithMultipleMeasurements(nrMeasurements);

		for (int c = 0; c < centroids.size(); ++c)
		{
			if (c == dataPoints[i].cluster) continue;

			const double dxo = dataPoints[i].x - centroids[c].x;
			const double dyo = dataPoints[i].y - centroids[c].y;
			const double d2o = dxo * dxo + dyo * dyo;

			if (d2o < d2 && d2 - d2o > 0.1)
			{
				std::cout << "Point " << i << " (" << dataPoints[i].x << ", " << dataPoints[i].y << ") is closer to centroid " << c << " (" << centroids[c].x << ", " << centroids[c].y << ") than to centroid " << dataPoints[i].cluster << " (" << centroids[dataPoints[i].cluster].x << ", " << centroids[dataPoints[i].cluster].y << ")" << std::endl;
				std::cout << "Euclidian: d2 for marked centroid: " << d2 << ", d2 for closer centroid: " << d2o << std::endl;
			}

			qMeansClustering.SetDataPoint2(centroids[c]);
			auto res2 = qMeansClustering.ExecuteWithMultipleMeasurements(nrMeasurements);
			if (res2[0] > res[0])
			{
				std::cout << "Point " << i << " (" << dataPoints[i].x << ", " << dataPoints[i].y << ") is closer to centroid " << c << " (" << centroids[c].x << ", " << centroids[c].y << ") than to centroid " << dataPoints[i].cluster << " (" << centroids[dataPoints[i].cluster].x << ", " << centroids[dataPoints[i].cluster].y << ") even with the quantum distance!" << std::endl;
				std::cout << "Measured count for marked centroid: " << res[0] << ", measured count for closer centroid: " << res2[0] << std::endl;
				return false;
			}
		}
	}

	return true;
}

bool QMLTests()
{
	return QMeansClustering2DTests();
}