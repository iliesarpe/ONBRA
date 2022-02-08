#include "graph.h"

#include <sstream>
#include <unordered_map>
#include <tuple>

// Removes timesteps (i.e. "compresses" the time) where there are no temporal edges
void removeBoringTimesteps(akt::TemporalEdgeSet& edges)
{
    int timestepsRemoved = 0, lastTimestep = 0;
    // Handle first edge separately
    if (!edges.empty())
        timestepsRemoved = edges.cbegin()->when;
    
    // Flatten the set into an array to be able to easily modify the when-values
    auto flat = std::vector<akt::TemporalEdge>(edges.cbegin(), edges.cend());
    
	// This is weird, why do we need this passage is just a waste of time and useless on very large networks
	// Why not to use better mapping in next steps?
	// Removing all empty timesteps (may be a problem in some applications)
    for (auto& te : flat) {
        te.when -= timestepsRemoved;
        if ((te.when - lastTimestep) > 1) {
            int newTime = lastTimestep + 1;
            timestepsRemoved += te.when - newTime;
            te.when = newTime;
        }
        lastTimestep = te.when;
    }
    
    // Re-blow into a tree
    edges = akt::TemporalEdgeSet(flat.cbegin(), flat.cend(), &akt::temporalEdgeLessTimewise);
}

namespace akt {
    bool temporalEdgeLessTimewise(const akt::TemporalEdge& lhs, const akt::TemporalEdge& rhs)
    {
        return  std::tie(lhs.when, lhs.from, lhs.to) < std::tie(rhs.when, rhs.from, rhs.to);
        //return (lhs.when != rhs.when) ? (lhs.when < rhs.when)
        //    : ((lhs.from != rhs.from) ? (lhs.from < rhs.from)
        //    : (lhs.to < rhs.to));
    }

    std::pair<Graph, std::vector<std::string>> readReduceGraph(std::istream& is, bool directed)
    {
        TemporalEdgeSet edgesRead(&temporalEdgeLessTimewise); // Passing the function pointer of the comparator
        std::unordered_map<std::string, int> nodeIds; // assuming at most max_int nodes
        std::vector<std::string> reverseNodeIds;
        int noNodes = 0;
        for (std::string line; std::getline(is, line); ) {
            auto iss = std::istringstream{ line };
            std::string from, to;
            int w;
            iss >> from >> to >> w;
			// REMAPPING NODES String -> int accoring to the first appearences
            if ((from.empty()) || (to.empty()))
                continue;
            if (nodeIds.count(from) < 1) {
                nodeIds[from] = noNodes++;
                reverseNodeIds.push_back(from);
            }
            if (nodeIds.count(to) < 1) {
                nodeIds[to] = noNodes++;
                reverseNodeIds.push_back(to);
            }
            int f = nodeIds[from];
            int t = nodeIds[to];

			//No self loops allowed
            if (f == t)
                continue;
            edgesRead.insert(TemporalEdge{ f, t, w });
			// IF graph is undirected then (0,1) means 0->1 and 1->0 and so on
            if (!directed)
                edgesRead.insert(TemporalEdge{ t, f, w });
        }
        removeBoringTimesteps(edgesRead);

        return { Graph(noNodes, edgesRead.crbegin()->when, edgesRead), reverseNodeIds };
    }
}
