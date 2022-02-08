#pragma once

#include <utility>
#include <vector>

#include "graph.h"
#include "adjList.h"

namespace akt {
	// Returns the shortest betweenness and shortest foremost betweenness (with strictness depending on the parameter)
	std::pair<std::vector<double>, std::vector<double>> shortestBetweenness(const akt::Graph& g, bool strict);

	// Returns the shortest betweenness approximated and shortest foremost betweenness (with strictness depending on the parameter) with adj
	std::pair<std::vector<double>, std::vector<double>> shortestBetweennessApproxAdj(const akt::GraphAdj& g, bool strict, int samples, int seed, double eta);

	std::pair<std::vector<double>, std::vector<double>> shortestBetweennessApproxAdjDelta(const akt::GraphAdj& g, bool strict, int samples, int seed, int delta, double eta);

} // end namespace akt
