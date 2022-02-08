#include "adjList.h"

#include <sstream>
#include <unordered_map>
#include <tuple>

/*
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
*/

namespace akt {
    bool temporalEdgeLessTimewise_64(const akt::TemporalEdge_64& lhs, const akt::TemporalEdge_64& rhs)
    {
        return  std::tie(lhs.when, lhs.from, lhs.to) < std::tie(rhs.when, rhs.from, rhs.to);
        //return (lhs.when != rhs.when) ? (lhs.when < rhs.when)
        //    : ((lhs.from != rhs.from) ? (lhs.from < rhs.from)
        //    : (lhs.to < rhs.to));
    }

    std::pair<GraphAdj, std::vector<std::string>> readReduceGraphToAdj(std::istream& is, bool directed)
    {
		// Using orderd set: therefore removing all (u,v,t) that are the same
        TemporalEdgeSet_64 edgesRead(&temporalEdgeLessTimewise_64); // Passing the function pointer of the comparator
        std::unordered_map<std::string, int> nodeIds; // assuming at most max_int nodes
        std::vector<std::string> reverseNodeIds;
        int noNodes = 0;
		long long int selfloops{0};
		int skipped{0};
		int inserted{0};
        for (std::string line; std::getline(is, line); ) 
		{
			//std::cout << line << std::endl;
            auto iss = std::istringstream{ line };
            std::string from, to;
            long long int w;
            iss >> from >> to >> w;
            if ((from.empty()) || (to.empty()))
			{
				skipped++;
                continue;
			}
			// Remapping source and dest to range in 0,...,n-1
            if (nodeIds.count(from) < 1) 
			{
                nodeIds[from] = noNodes++;
                reverseNodeIds.push_back(from);
            }
            if (nodeIds.count(to) < 1) 
			{
                nodeIds[to] = noNodes++;
                reverseNodeIds.push_back(to);
            }
			// REMAPPING NODES String -> int accoring to the first appearences
            int f = nodeIds[from];
            int t = nodeIds[to];

			//No self loops allowed
            if (f == t)
			{
				selfloops++;
                continue;
			}
            auto ret = edgesRead.insert(TemporalEdge_64{ f, t, w });
			if(ret.second)
				inserted++;
			// IF graph is undirected then (0,1) means 0->1 and 1->0 and so on
            if (!directed)
                edgesRead.insert(TemporalEdge_64{ t, f, w });
		}
		std::cout << "Self-loops removed: "<< selfloops << std::endl;
		std::cout << "Skipped: "<< skipped << std::endl;
		std::cout << "Edges inserted: "<< inserted << std::endl;
        return { GraphAdj(noNodes, edgesRead.crbegin()->when, std::cref(edgesRead)), reverseNodeIds };
    }
}
