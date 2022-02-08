#include "algorithms.h"

#include <chrono>
#include <algorithm>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <random>
#include <chrono>
#include <boost/functional/hash.hpp>

bool temporalEdgeGreaterTimewise(const akt::TemporalEdge& lhs, const akt::TemporalEdge& rhs)
{
    return (lhs.when != rhs.when) ? (lhs.when > rhs.when)
        : ((lhs.from != rhs.from) ? (lhs.from < rhs.from)
            : (lhs.to < rhs.to));
}

// Stores the results of a bfs into a temporal graph: the shortest distance to it and the earliest possible arrival time there
struct VertexDistInfo
{
    // The shortest distance to the node
    int dist;
    // The earliest possible arrival time at the node
    int foremostTime;
};

// Stores data about a vertex appearance
struct VertexAppearance
{
    // Vertex id
    int v;
    // Timestamp
    long long int time;
};
 
bool operator==(const VertexAppearance& lhs, const VertexAppearance& rhs)
{
    return (lhs.v == rhs.v) && (lhs.time == rhs.time);
}

// Stores data about a vertex appearance
struct VertexAppearanceSuccint
{
    // Vertex id
    int v;
    // Timestamp
    long long int time;
    // Index of this vertex appeareance
    long long int idx;
};


bool operator==(const VertexAppearanceSuccint& lhs, const VertexAppearanceSuccint& rhs)
{
    return (lhs.v == rhs.v) && (lhs.time == rhs.time);
}

// Simple hashing function for VertexAppearance (so that we can use it with std::unordered_set)
// Adapted from https://en.cppreference.com/w/cpp/utility/hash
// Better hash functions using Knuth's observation
namespace std
{
    template<> struct hash<VertexAppearance>
    {
        std::size_t operator()(const VertexAppearance& va) const noexcept
        {
            auto h1 = std::hash<long long int>{}(va.v);
            auto h2 = std::hash<long long int>{}(va.time);
            return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };

	template<> struct hash<VertexAppearanceSuccint>
    {
        std::size_t operator()(const VertexAppearanceSuccint& va) const noexcept
        {
			std::size_t seed = 0;
			boost::hash_combine(seed, va.v * 2654435761);
			boost::hash_combine(seed, va.time * 2654435761);
			boost::hash_combine(seed, va.idx * 2654435761);
			return seed;
            //auto h1 = std::hash<long long int>{}(va.v);
            //auto h2 = std::hash<long long int>{}(va.time);
            //return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };
}

// Struct for storing the helper arrays for each iteration of the outer loop of the algorithm for betweenness centrality
struct ShortestBetweennessData
{
	// TODO: Maybe rephrase all of this in a memory efficient way?
    ShortestBetweennessData(int n, int T)
        : deltaDots{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) },
          sigmas{ std::vector<std::vector<int>>(n, std::vector<int>(T, 0)) },
          preds{ std::vector<std::vector<std::unordered_set<VertexAppearance>>>(n, std::vector<std::unordered_set<VertexAppearance>>(T)) },
          dists{ std::vector<std::vector<int>>(n, std::vector<int>(T, -1)) },
          totalDists{ std::vector<int>(n, -1) },
          totalSigmas{ std::vector<int>(n, 0) },
          stack{ std::vector<VertexAppearance>() },
          foremostTimes{ std::vector<int>(n, -1) },
          deltaForemosts{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) }
    { }
          
    ShortestBetweennessData(const akt::Graph& g)
        : ShortestBetweennessData(g.N(), g.T())
    { }

    // The \delta_{s\cdot} array on each iteration
    std::vector<std::vector<double>> deltaDots;

    // Numbers of shortest paths from a source to each vertex appearance
    std::vector<std::vector<int>> sigmas;
	
    // Sets of predecessors for each vertex appearance
    std::vector<std::vector<std::unordered_set<VertexAppearance>>> preds;

    // Distances to each vertex appearance
    std::vector<std::vector<int>> dists;

    // Distances to each vertex (*not* vertex __appearance__)
    std::vector<int> totalDists;

    // Number of shortests paths from a source to each vertex (*not* vertex __appearance__)
    std::vector<int> totalSigmas;

    // Stack of vertex apperances in order of their discovery with the bfs
    std::vector<VertexAppearance> stack;

    // Additional arrays for handling shortest foremost computation
    // Earliest arrival times for each vertex
    std::vector<int> foremostTimes;
    // Delta array for foremost paths
    std::vector<std::vector<double>> deltaForemosts;
};

// Struct for storing the helper arrays for each iteration of the outer loop of the algorithm for betweenness centrality - Succint
struct ShortestBetweennessDataSuccint
{
	// TODO: Maybe rephrase all of this in a memory efficient way?
    ShortestBetweennessDataSuccint(int n, int T, const std::vector<long long int>& outtimes, const std::vector<std::vector<akt::GraphAdj::AppearanceNeighbourhood>>& adjacencyList)
        : //deltaDots{ std::vector<std::vector<double>>(n, std::vector<double>()) },
          sigmas{ std::vector<std::vector<std::pair<long long int, long long int> >>(n, std::vector<std::pair<long long int,long long  int> >()) },
          //sigmas{ std::vector<std::vector<std::pair<long long int, double> >>(n, std::vector<std::pair<long long int,double> >()) },
          preds{ std::vector<std::vector<std::pair<long long int, std::unordered_set<VertexAppearanceSuccint>> >>(n, std::vector<std::pair<long long int, std::unordered_set<VertexAppearanceSuccint>> >()) },
          //preds{ std::vector<std::vector<std::pair<long long int, std::unordered_set<VertexAppearance>> >>(n, std::vector<std::pair<long long int, std::unordered_set<VertexAppearance>> >()) },
          dists{ std::vector<std::vector<std::pair<long long int, long long int> >>(n, std::vector<std::pair<long long int,long long  int> >()) },
          totalDists{ std::vector<long long int>(n, -1) },
          totalSigmas{ std::vector<long long int>(n, 0) },
		  stack{ std::vector<VertexAppearanceSuccint>() },
          foremostTimes{ std::vector< long long int>(n, -1) },
          //stack{ std::vector<VertexAppearance>() },
          sigmas_to_z{ std::vector<std::vector<std::pair<long long int, long long int> >>(n, std::vector<std::pair<long long int,long long  int> >()) },
          //sigmas_to_z{ std::vector<std::vector<std::pair<long long int, double> >>(n, std::vector<std::pair<long long int,double> >()) },
          inQueue{ std::vector<std::vector<std::pair<long long int, bool> >>(n, std::vector<std::pair<long long int, bool> >()) },
          inQueueUsages{ std::vector<std::vector<std::pair<long long int, long long int> >>(n, std::vector<std::pair<long long int, long long int> >()) },
          visitedNodes{ std::vector<bool>(n, true) }
          //deltaForemosts{ std::vector<std::vector<double>>(n, std::vector<double>(T, 0.0)) }
    { 
		size_t idx=0;
		// In case of source
		for(size_t i{0}; i < static_cast<size_t>(n); i++)
		{
			sigmas[i].push_back(std::make_pair(0, 0));
			//sigmas[i].push_back(std::make_pair(0, 0.));
			dists[i].push_back(std::make_pair(0, -1));
			preds[i].push_back(std::make_pair(0, std::unordered_set<VertexAppearanceSuccint>()));
			//preds[i].push_back(std::make_pair(0, std::unordered_set<VertexAppearance>()));
			sigmas_to_z[i].push_back(std::make_pair(0, 0));
			//sigmas_to_z[i].push_back(std::make_pair(0, 0.));
			inQueue[i].push_back(std::make_pair(0, false));
			inQueueUsages[i].push_back(std::make_pair(0, 0));
		}
		for(auto& entry: adjacencyList)	// For each node
		{
			for(auto& out : entry) // For each outgoing timestamp
			{
				for(auto& neigh : out.neighbours) // for each neighbour at the given timestamp
				{
					//if((timestamps.count(out.timestep)==0) || (sigmas[neigh].size()==0))
					//{
				//		sigmas[neigh].push_back(std::make_pair(out.timestep, 0));
				//		dists[neigh].push_back(std::make_pair(out.timestep, -1));
				//		preds[neigh].push_back(std::make_pair(out.timestep, std::unordered_set<VertexAppearance>()));
				//	}
					//else // Timestamp already seen
					//{
						auto pos = std::lower_bound(sigmas[neigh].begin(), sigmas[neigh].end(), out.timestep, [](auto& pair, auto tim){return (pair.first < tim);}) - sigmas[neigh].begin();
						//std::cout << "pos is: " << pos << " size: " << sigmas[neigh].size() << " tim: " << out.timestep << std::endl;
						if(pos==static_cast<long long int>(sigmas[neigh].size())||(sigmas[neigh].size()==1))
						{
							sigmas[neigh].push_back(std::make_pair(out.timestep, 0));
							dists[neigh].push_back(std::make_pair(out.timestep, -1));
							preds[neigh].push_back(std::make_pair(out.timestep, std::unordered_set<VertexAppearanceSuccint>()));
							//preds[neigh].push_back(std::make_pair(out.timestep, std::unordered_set<VertexAppearance>()));
							sigmas_to_z[neigh].push_back(std::make_pair(out.timestep, 0));
							inQueue[neigh].push_back(std::make_pair(out.timestep, false));
							inQueueUsages[neigh].push_back(std::make_pair(out.timestep, 0));
						}
						else if(sigmas[neigh][pos].first > out.timestep)
						{
							sigmas[neigh].insert(sigmas[neigh].begin()+pos, std::make_pair(out.timestep, 0));
							dists[neigh].insert(dists[neigh].begin()+pos, std::make_pair(out.timestep, -1));
							preds[neigh].insert(preds[neigh].begin()+pos, std::make_pair(out.timestep, std::unordered_set<VertexAppearanceSuccint>()));
							//preds[neigh].insert(preds[neigh].begin()+pos, std::make_pair(out.timestep, std::unordered_set<VertexAppearance>()));
							sigmas_to_z[neigh].insert(sigmas_to_z[neigh].begin()+pos, std::make_pair(out.timestep, 0));
							inQueue[neigh].insert(inQueue[neigh].begin()+pos,std::make_pair(out.timestep, false));
							inQueueUsages[neigh].insert(inQueueUsages[neigh].begin()+pos,std::make_pair(out.timestep, 0));
						}
						/*
						if(pos == sigmas[neigh].size())
						{
						}
						else
						{ || sigmas[neigh][pos-1].first > out.timestep) || (sigmas[neigh][pos-1].first > out.timestep)
						}
						*/
					//}
				}
			}
			idx++;
		}
		std::cout << " Ended " << std::endl;
	}
          
    ShortestBetweennessDataSuccint(const akt::GraphAdj& g)
        : ShortestBetweennessDataSuccint(g.N(), g.T(), g.getOutTimestamps(), g.adjacencyList()) // g.Timeline();
    { }

    // The \delta_{s\cdot} array on each iteration
    std::vector<std::vector<double>> deltaDots;

    // Numbers of shortest paths from a source to each vertex appearance
    std::vector<std::vector<std::pair<long long int,long long  int> >> sigmas;
    //std::vector<std::vector<std::pair<long long int,double> >> sigmas;
	
    // Sets of predecessors for each vertex appearance
    std::vector<std::vector<std::pair<long long int, std::unordered_set<VertexAppearanceSuccint>>>> preds;
    //std::vector<std::vector<std::pair<long long int, std::unordered_set<VertexAppearance>>>> preds;

    // Distances to each vertex appearance
    std::vector<std::vector<std::pair<long long int,long long  int> >> dists;

    // Distances to each vertex (*not* vertex __appearance__)
    std::vector<long long int> totalDists;

    // Number of shortests paths from a source to each vertex (*not* vertex __appearance__)
    std::vector<long long int> totalSigmas;

    // Stack of vertex apperances in order of their discovery with the bfs
    std::vector<VertexAppearanceSuccint> stack;
    //std::vector<VertexAppearance> stack;

    // Additional arrays for handling shortest foremost computation
    // Earliest arrival times for each vertex
    std::vector<long long int> foremostTimes;
    // Delta array for foremost paths
    std::vector<std::vector<double>> deltaForemosts;

	// Sigmas to z used in computing backwards betweenness
	std::vector<std::vector<std::pair<long long int, long long int>>> sigmas_to_z;
	//std::vector<std::vector<std::pair<long long int, double>>> sigmas_to_z;

	// Is the vertex appeareance already visited when computing backwards centralities?
	std::vector<std::vector<std::pair<long long int, bool>>> inQueue;

	//How many times Is the vertex appeareance already visited when computing backwards centralities? WARN: used only in the delta-variant
	std::vector<std::vector<std::pair<long long int, long long int>>> inQueueUsages;

	// Vector that keeps track for each nde if it was explored during the BFS
    std::vector<bool> visitedNodes;

};

void reinitializeHelperStructSuccint(const akt::GraphAdj& g, int s, ShortestBetweennessDataSuccint& sbd, int prevsource)
{
	/*
	if(prevsource > 0)
	{
		//std::fill(sbd.sigmas[prevsource].begin(), sbd.sigmas[prevsource].end(), 0);
		//std::fill(sbd.dists[prevsource].begin(), sbd.dists[prevsource].end(), -1);
		//std::fill(sbd.deltaDots[prevsource].begin(), sbd.deltaDots[prevsource].end(), 0.0);
		//std::fill(sbd.deltaForemosts[prevsource].begin(), sbd.deltaForemosts[prevsource].end(), 0.0);
		std::for_each(sdb.sigmas[prevsource].begin(), sbd.sigmas[prevsource].end(), [](auto& entry) {}
	}
    // Reinitialize the appearance-based arrays (at the only relevant times)
    for (auto it = g.edges_cbegin(); it != g.edges_cend(); ++it) {
        auto te = *it;
        sbd.deltaDots[te.to][te.when] = 0.0;
        sbd.deltaForemosts[te.to][te.when] = 0.0;
        sbd.sigmas[te.to][te.when] = 0;
        sbd.preds[te.to][te.when].clear();
        sbd.dists[te.to][te.when] = -1;
    }
	*/
	//long long int prev{0};
	for(size_t i{0}; i < static_cast<size_t>(g.N()); i++)
	{
		// Reset only if visited
		if(sbd.visitedNodes[i])
		{
			//for(long long int j{0}; j < g.getOutTimestamps()[i]; j++)
			for(size_t j{0}; j < sbd.sigmas[i].size(); j++)
			{
				//if(j>0)
				//{
					//std::cout << "noode " << i << " indeg time pos: " << j << " tim: " << sbd.sigmas[i][j].first << std::endl;
					//if(prev>sbd.sigmas[i][j].first)
					//	std::cout << "ERROR, vecotr not sorted!!" << i << " " << j << " " << prev << " " << sbd.sigmas[i][j].first << '\n';
				//}
				//prev=sbd.sigmas[i][j].first;
				sbd.sigmas[i][j].second = 0;
				//sbd.sigmas[i][j].second = 0.;
				sbd.preds[i][j].second.clear();
				sbd.dists[i][j].second = -1;
				sbd.sigmas_to_z[i][j].second = 0;
				//sbd.sigmas_to_z[i][j].second = 0.;
				sbd.inQueue[i][j].second = false;
				sbd.inQueueUsages[i][j].second = 0;
			}
		}
	}
    // Reinitialize the vertex-based arrays
    std::transform(sbd.totalDists.cbegin(), sbd.totalDists.cend(), sbd.totalDists.begin(), [](auto i) { return -1; });
    std::transform(sbd.totalSigmas.cbegin(), sbd.totalSigmas.cend(), sbd.totalSigmas.begin(), [](auto i) { return 0; });
		//std::for_each(sbd.totalSigmas.begin(), sbd.totalSigmas.end(), [](auto& el) {el=0;});
    std::transform(sbd.foremostTimes.cbegin(), sbd.foremostTimes.cend(), sbd.foremostTimes.begin(), [](auto i) { return -1; });
    std::transform(sbd.visitedNodes.cbegin(), sbd.visitedNodes.cend(), sbd.visitedNodes.begin(), [](auto i) { return false; });
    // Reinitialize the stack (though it should be empty already, so it's just a sanity operation)
    sbd.stack.clear();
    // Initialize the elements involving the source to appropriate values (if different from the defaults)
	sbd.visitedNodes[s] = true;
    sbd.sigmas[s][0].second = 1;
    //sbd.sigmas[s][0].second = 1.;
    sbd.dists[s][0].second = 0;
    sbd.totalDists[s] = 0;
    sbd.totalSigmas[s] = 1;
    sbd.foremostTimes[s] = 0;
}

// Struct for storing the helper arrays for the prefix-foremost computation algorithm
struct PrefixBetweennessData
{
    PrefixBetweennessData(int n)
        : deltaDots(n, 1.0), sigmas(n, 0), preds(n), foremostTimes(n, -1)
    { }

    PrefixBetweennessData(const akt::Graph& g)
        : PrefixBetweennessData(g.N())
    { }

    // The \delta_{s\cdot} values for the foremost apperances
    std::vector<double> deltaDots;
    // Numbers of prefix-foremost paths from the source to each vertex
    std::vector<int> sigmas;
    // Sets of predecessors for each vertex
    std::vector<std::unordered_set<int>> preds;
    // foremost time to each vertex
    std::vector<int> foremostTimes;
    // Stack of nodes in order of their discovery with the "foremost-based" search
    std::vector<int> stack;
};

// (Re-) Initializes all the members in sbd for the next iteration of the outermost iteration of the shortest betwenness algorithm
void reinitializeHelperStruct(const akt::Graph& g, int s, ShortestBetweennessData& sbd, int prevsource)
{
	if(prevsource > 0)
	{
		/*
    	for(auto& value : sbd.sigmas[prevsource])
			value = 0;
    	for(auto& value : sbd.dists[prevsource])
			value = -1;
    	for(auto& value : sbd.deltaDots[prevsource])
        	value = 0.0;
    	for(auto& value : sbd.deltaForemosts[prevsource])
        	value = 0.0;
		*/
		std::fill(sbd.sigmas[prevsource].begin(), sbd.sigmas[prevsource].end(), 0);
		std::fill(sbd.dists[prevsource].begin(), sbd.dists[prevsource].end(), -1);
		std::fill(sbd.deltaDots[prevsource].begin(), sbd.deltaDots[prevsource].end(), 0.0);
		std::fill(sbd.deltaForemosts[prevsource].begin(), sbd.deltaForemosts[prevsource].end(), 0.0);
	}
    // Reinitialize the appearance-based arrays (at the only relevant times)
    for (auto it = g.edges_cbegin(); it != g.edges_cend(); ++it) {
        auto te = *it;
        sbd.deltaDots[te.to][te.when] = 0.0;
        sbd.deltaForemosts[te.to][te.when] = 0.0;
        sbd.sigmas[te.to][te.when] = 0;
        sbd.preds[te.to][te.when].clear();
        sbd.dists[te.to][te.when] = -1;
    }
    // Reinitialize the vertex-based arrays
    std::transform(sbd.totalDists.cbegin(), sbd.totalDists.cend(), sbd.totalDists.begin(), [](auto i) { return -1; });
    std::transform(sbd.totalSigmas.cbegin(), sbd.totalSigmas.cend(), sbd.totalSigmas.begin(), [](auto i) { return 0; });
		//std::for_each(sbd.totalSigmas.begin(), sbd.totalSigmas.end(), [](auto& el) {el=0;});
    std::transform(sbd.foremostTimes.cbegin(), sbd.foremostTimes.cend(), sbd.foremostTimes.begin(), [](auto i) { return -1; });
    // Reinitialize the stack (though it should be empty already, so it's just a sanity operation)
    sbd.stack.clear();
    // Initialize the elements involving the source to appropriate values (if different from the defaults)
    sbd.sigmas[s][0] = 1;
    sbd.dists[s][0] = 0;
    sbd.totalDists[s] = 0;
    sbd.totalSigmas[s] = 1;
    sbd.foremostTimes[s] = 0;
}

// (Re-) Initializes all the members in sbd for the next iteration of the outermost iteration of the shortest betwenness algorithm
void reinitializeHelperStructSampling(const akt::Graph& g, int s, ShortestBetweennessData& sbd, int prevS)
{
    // Reinitialize the appearance-based arrays (at the only relevant times)
	/*
    for (auto it = g.edges_cbegin(); it != g.edges_cend(); ++it) {
        auto te = *it;
        sbd.deltaDots[te.to][te.when] = 0.0;
        sbd.deltaForemosts[te.to][te.when] = 0.0;
        sbd.sigmas[te.to][te.when] = 0;
        sbd.preds[te.to][te.when].clear();
        sbd.dists[te.to][te.when] = -1;
    }
	*/
	if(prevS > 0)
	{
    	for(auto& value : sbd.sigmas[prevS])
			value = 0;
    	for(auto& value : sbd.dists[prevS])
			value = -1;
	}
			
	for(auto& appeareance : sbd.stack)
	{
        //sbd.deltaDots[appeareance.v][appeareance.time] = 0.0;
        //sbd.deltaForemosts[appeareance.v][appeareance.time] = 0.0;
        sbd.sigmas[appeareance.v][appeareance.time] = 0;
        sbd.preds[appeareance.v][appeareance.time].clear();
        sbd.dists[appeareance.v][appeareance.time] = -1;
	}
    // Reinitialize the vertex-based arrays
    std::transform(sbd.totalDists.cbegin(), sbd.totalDists.cend(), sbd.totalDists.begin(), [](auto i) { return -1; });
    std::transform(sbd.totalSigmas.cbegin(), sbd.totalSigmas.cend(), sbd.totalSigmas.begin(), [](auto i) { return 0; });
		//std::for_each(sbd.totalSigmas.begin(), sbd.totalSigmas.end(), [](auto& el) {el=0;});
    //std::transform(sbd.foremostTimes.cbegin(), sbd.foremostTimes.cend(), sbd.foremostTimes.begin(), [](auto i) { return -1; });
    // Reinitialize the stack (though it should be empty already, so it's just a sanity operation)
    sbd.stack.clear();
    // Initialize the elements involving the source to appropriate values (if different from the defaults)
    sbd.sigmas[s][0] = 1;
    sbd.dists[s][0] = 0;
    sbd.totalDists[s] = 0;
    sbd.totalSigmas[s] = 1;
    //sbd.foremostTimes[s] = 0;
}

// Computes all the distances to vertices, vertex appearances and the counts of corresponding shortest paths
void shortestComputeDistancesSigmas(const akt::Graph& g, bool strict, int s, ShortestBetweennessData& sbd)
{
    std::queue<VertexAppearance> q;
    q.push(VertexAppearance{ s, 0 });
    // BFS to find distances and sigma values
    while (!q.empty()) {
        auto cur = q.front();
        q.pop();
        // Go over all neighbours of the current vertex appearance
        for (int t = cur.time + strict; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep) {
            for (auto w : g.adjacencyList()[cur.v][t].neighbours) {
                if (sbd.dists[w][t] < 0) {
                    sbd.dists[w][t] = sbd.dists[cur.v][cur.time] + 1;
					// We may not need this
                    if (sbd.totalDists[w] < 0)
                        sbd.totalDists[w] = sbd.dists[w][t];
                    q.push(VertexAppearance{ w, t });
                    sbd.stack.push_back(VertexAppearance{ w, t });
                }
                if (sbd.dists[w][t] == sbd.dists[cur.v][cur.time] + 1) {
                    sbd.sigmas[w][t] += sbd.sigmas[cur.v][cur.time];
                    sbd.preds[w][t].insert(cur);
                    if (sbd.totalDists[w] == sbd.dists[cur.v][cur.time] + 1)
                        sbd.totalSigmas[w] += sbd.sigmas[cur.v][cur.time];
                }
                if ((sbd.foremostTimes[w] < 0) || (t < sbd.foremostTimes[w])) {
                    sbd.foremostTimes[w] = t;
                }
            }
        }
    }
}

// Computes all the distances to vertices, vertex appearances and the counts of corresponding shortest paths
bool shortestComputeDistancesSigmasTruncated(const akt::Graph& g, bool strict, int s, int z, ShortestBetweennessData& sbd)
{
    std::queue<VertexAppearance> q;
	int mindisttoz{-1};
    q.push(VertexAppearance{ s, 0 });
    // BFS to find distances and sigma values
    while (!q.empty()) {
        auto cur = q.front();
        q.pop();
		if(((mindisttoz > 0) && sbd.dists[cur.v][cur.time] < mindisttoz) || mindisttoz == -1)
		{
        // Go over all neighbours of the current vertex appearance
        for (int t = cur.time + strict; (t >= 0) && (t <= g.maximalTimestep()); t = g.adjacencyList()[cur.v][t].nextTimestep) {
            for (auto w : g.adjacencyList()[cur.v][t].neighbours) {
                if (sbd.dists[w][t] < 0) {
                    sbd.dists[w][t] = sbd.dists[cur.v][cur.time] + 1;
					// First time seeing node w
                    if (sbd.totalDists[w] < 0)
					{
                        sbd.totalDists[w] = sbd.dists[w][t];
						if(w == z)
							mindisttoz = sbd.dists[w][t];
					}
					//if reached z do not push push the current vertex appeareance as candidate
					if(w != z)
					{
						if((mindisttoz < 0) || ((mindisttoz > 0) && (sbd.dists[w][t] < mindisttoz)))
                    		q.push(VertexAppearance{ w, t });
					}
                    sbd.stack.push_back(VertexAppearance{ w, t });
                }
                if (sbd.dists[w][t] == sbd.dists[cur.v][cur.time] + 1) {
                    sbd.sigmas[w][t] += sbd.sigmas[cur.v][cur.time];
                    sbd.preds[w][t].insert(cur);
					//This may not be needed
                    if (sbd.totalDists[w] == sbd.dists[cur.v][cur.time] + 1)
                        sbd.totalSigmas[w] += sbd.sigmas[cur.v][cur.time];
                }
            }
        }
		}
    }
	return mindisttoz > 0;
}

bool shortestComputeDistancesSigmasTruncatedSuccint(const akt::GraphAdj& g, bool strict, int s, int z, ShortestBetweennessDataSuccint& sbd)
{
    std::queue<VertexAppearanceSuccint> q;
	int mindisttoz{-1};
    q.push(VertexAppearanceSuccint{ s, 0, 0 });
	const std::vector<long long int>& outtim = g.getOutTimestamps();
	auto adj = g.adjacencyList();
    // BFS to find distances and sigma values
    while (!q.empty()) {
        auto cur = q.front();
        q.pop();
		auto posAdj = std::lower_bound(adj[cur.v].begin(), adj[cur.v].end(), cur.time+strict, [](auto& entry, long long int tim){ return (entry.timestep < tim) ;}) - adj[cur.v].begin();
		auto posDist = cur.idx;
		if(((mindisttoz > 0) && sbd.dists[cur.v][posDist].second < mindisttoz) || (mindisttoz == -1))
		{
			// Go over all neighbours of the current vertex appearance
			for(auto index{posAdj}; index < outtim[cur.v]; index++)
			//for(auto index{posAdj}; index < adj[cur.v].size(); index++)
			{
				long long int t = adj[cur.v][index].timestep;
				// check if t is within at most delta time, otherwise prune this path
				for (auto w : adj[cur.v][index].neighbours) {
					sbd.visitedNodes[w] = true;
					auto post = std::lower_bound(sbd.dists[w].begin(), sbd.dists[w].end(), t, [](auto& pair, long long int tim) {return (pair.first < tim);}) - sbd.dists[w].begin();
					if (sbd.dists[w][post].second < 0) {

						sbd.dists[w][post].second = sbd.dists[cur.v][posDist].second + 1;
						// First time seeing node w
						if (sbd.totalDists[w] < 0)
						{
							sbd.totalDists[w] = sbd.dists[w][post].second;
							if(w == z)
							{
								mindisttoz = sbd.dists[w][post].second;
							}
						}
						//if reached z do not push push the current vertex appeareance as candidate
						if(w != z)
						{
							if((mindisttoz < 0) || ((mindisttoz > 0) && (sbd.dists[w][post].second < mindisttoz)))
							{
								q.push(VertexAppearanceSuccint{ w, t, post });
							}
						}
						sbd.stack.push_back(VertexAppearanceSuccint{ w, t, post });
					}
					if (sbd.dists[w][post].second == sbd.dists[cur.v][posDist].second + 1) {
						sbd.sigmas[w][post].second += sbd.sigmas[cur.v][posDist].second;
						sbd.preds[w][post].second.insert(cur);
						if (sbd.totalDists[w] == sbd.dists[cur.v][posDist].second + 1)
							sbd.totalSigmas[w] += sbd.sigmas[cur.v][posDist].second;
					}
				}
			}
		}
    }

	return mindisttoz > 0;
}

bool shortestComputeDistancesSigmasTruncatedSuccintDelta(const akt::GraphAdj& g, bool strict, int s, int z, ShortestBetweennessDataSuccint& sbd, int delta)
{
	//TODO: Still to truncate computatation
    std::queue<VertexAppearanceSuccint> q;
	int mindisttoz{-1};
    q.push(VertexAppearanceSuccint{ s, 0, 0 });
	const std::vector<long long int>& outtim = g.getOutTimestamps();
	auto adj = g.adjacencyList();
    // BFS to find distances and sigma values
    while (!q.empty()) {
        auto cur = q.front();
        q.pop();
		auto posAdj = std::lower_bound(adj[cur.v].begin(), adj[cur.v].end(), cur.time+strict, [](auto& entry, long long int tim){ return (entry.timestep < tim) ;}) - adj[cur.v].begin();
		auto posDist = cur.idx;
		if(((mindisttoz > 0) && sbd.dists[cur.v][posDist].second < mindisttoz) || (mindisttoz == -1))
		{
			// Go over all neighbours of the current vertex appearance
			for(auto index{posAdj}; index < outtim[cur.v]; index++)
			{
				long long int t = adj[cur.v][index].timestep;
				// RESTLESS-VARIANT
				if((cur.time == 0) || (cur.time+delta >= t))
				{
				// check if t is within at most delta time, otherwise prune this path
				for (auto w : adj[cur.v][index].neighbours) {
					sbd.visitedNodes[w] = true;
					auto post = std::lower_bound(sbd.dists[w].begin(), sbd.dists[w].end(), t, [](auto& pair, long long int tim) {return (pair.first < tim);}) - sbd.dists[w].begin();
					//std::cout <<  " " << sbd.dists[w][post].first << " " <<  sbd.dists[w][post].second << std::endl;
					if (sbd.dists[w][post].second < 0) {

						sbd.dists[w][post].second = sbd.dists[cur.v][posDist].second + 1;
						// First time seeing node w
						if (sbd.totalDists[w] < 0)
						{
							sbd.totalDists[w] = sbd.dists[w][post].second;
							if(w == z)
							{
								mindisttoz = sbd.dists[w][post].second;
							}
						}
						//if reached z do not push push the current vertex appeareance as candidate
						if(w != z)
						{
							if((mindisttoz < 0) || ((mindisttoz > 0) && (sbd.dists[w][post].second < mindisttoz)))
							{
								q.push(VertexAppearanceSuccint{ w, t, post });
							}
						}
						sbd.stack.push_back(VertexAppearanceSuccint{ w, t, post });
					}
					if (sbd.dists[w][post].second == sbd.dists[cur.v][posDist].second + 1) {
						sbd.sigmas[w][post].second += sbd.sigmas[cur.v][posDist].second;
						sbd.preds[w][post].second.insert(cur);
						if (sbd.totalDists[w] == sbd.dists[cur.v][posDist].second + 1)
							sbd.totalSigmas[w] += sbd.sigmas[cur.v][posDist].second;
					}
				}
				}
			}
		}
    }

	return mindisttoz > 0;
}

void updateEstimates(int s, int z, ShortestBetweennessData& sbd, std::vector<double>& betweenness, const akt::Graph& g)
{
	std::vector<std::vector<double>> sigmas_to_z(betweenness.size(), std::vector<double>(g.T(), 0.));
	std::vector<std::vector<bool>> inQueue(betweenness.size(), std::vector<bool>(g.T(), false));
	std::vector<double> currentUpdateBetweenness(betweenness.size(), 0.);
    // Compute the delta array values (bottom-up) and update the betweenness values
    // Subtract the sum of the s-th row of the connectivity matrix to betweenness of s (follows from our formula)
	long long int numPathssz{sbd.totalSigmas[z]};
	std::queue<VertexAppearance> reverseVisit;
	// For each time that reaches the destination
	for(int time{static_cast<int>(sbd.preds[z].size())-1}; time >= 0; time--)
	{
		//If there exists a path that reaches the vertex at that time
		if(sbd.preds[z][time].size() > 0)
		{
			//Span each predecessor reaching the node at that time
			for (auto pred : sbd.preds[z][time]) 
			{
				sigmas_to_z[pred.v][pred.time] += 1;
				//Each predecessor has to be "visited backwards" only once
				if(inQueue[pred.v][pred.time] == false)
				{
					reverseVisit.push(VertexAppearance{ pred.v, pred.time });
					inQueue[pred.v][pred.time] = true;
				}
			}
		}
	}

	// While not computing all the betweenness values
    while (!reverseVisit.empty()) {
        auto cur = reverseVisit.front();
        reverseVisit.pop();
		if(cur.v == s)
		{
		}
		else
		{
			double pathstoz = sigmas_to_z[cur.v][cur.time];
			currentUpdateBetweenness[cur.v] += sbd.sigmas[cur.v][cur.time] * pathstoz;

        	for (auto pred : sbd.preds[cur.v][cur.time]) 
			{
				sigmas_to_z[pred.v][pred.time] += pathstoz;
				if(inQueue[pred.v][pred.time] == false)
				{
					reverseVisit.push(VertexAppearance{ pred.v, pred.time });
					inQueue[pred.v][pred.time] = true;
				}
			}
		}
    }

	for(size_t j{0}; j < currentUpdateBetweenness.size(); j++)
	{
		if(currentUpdateBetweenness[j] > 0)
			betweenness[j] += currentUpdateBetweenness[j]/numPathssz;

	}
	
}

void updateEstimatesSuccint(int s, int z, ShortestBetweennessDataSuccint& sbd, std::vector<double>& betweenness, const akt::GraphAdj& g)
{
	std::vector<double> currentUpdateBetweenness(betweenness.size(), 0.);
    // Compute the delta array values (bottom-up) and update the betweenness values
    // Subtract the sum of the s-th row of the connectivity matrix to betweenness of s (follows from our formula)
	long long int numPathssz{sbd.totalSigmas[z]};
	std::queue<VertexAppearanceSuccint> reverseVisit;
	// For each time that reaches the destination
	for(long long int idx_time{static_cast<long long int>(sbd.preds[z].size())-1}; idx_time >= 0; idx_time--)
	{
		//If there exists a path that reaches the vertex at that time
		if(sbd.preds[z][idx_time].second.size() > 0)
		{
			for (auto pred : sbd.preds[z][idx_time].second) 
			{
				auto cur_idx_time = pred.idx;
				sbd.sigmas_to_z[pred.v][cur_idx_time].second += 1;
				//Each predecessor has to be "visited backwards" only once
				if(sbd.inQueue[pred.v][cur_idx_time].second == false)
				{
					reverseVisit.push(VertexAppearanceSuccint{ pred.v, pred.time, pred.idx });
					sbd.inQueue[pred.v][cur_idx_time].second = true;
				}
			}
		}
	}

	bool f{true};
	// While not computing all the betweenness values
    while (!reverseVisit.empty()) {
        auto cur = reverseVisit.front();
        reverseVisit.pop();
		
		if(cur.v == s)
		{
		}
		else
		{
			auto cur_idx_time = cur.idx;
			auto pathstoz = sbd.sigmas_to_z[cur.v][cur_idx_time].second;
			

			currentUpdateBetweenness[cur.v] += sbd.sigmas[cur.v][cur_idx_time].second * pathstoz;

        	for (auto pred : sbd.preds[cur.v][cur_idx_time].second) 
			{
				auto pred_idx_time = pred.idx;
				sbd.sigmas_to_z[pred.v][pred_idx_time].second += pathstoz;
				if(sbd.inQueue[pred.v][pred_idx_time].second == false)
				{
					reverseVisit.push(VertexAppearanceSuccint{ pred.v, pred.time, pred.idx });
					sbd.inQueue[pred.v][pred_idx_time].second = true;
				}
			}
		}
    }

	for(size_t j{0}; j < currentUpdateBetweenness.size(); j++)
	{
		if(currentUpdateBetweenness[j] > 0)
		{
			betweenness[j] += currentUpdateBetweenness[j]/static_cast<double>(numPathssz);
		}

	}
}

void updateEstimatesSuccintDelta(int s, int z, ShortestBetweennessDataSuccint& sbd, std::vector<double>& betweenness, const akt::GraphAdj& g)
{
	std::vector<double> currentUpdateBetweenness(betweenness.size(), 0.);
    // Compute the delta array values (bottom-up) and update the betweenness values
    // Subtract the sum of the s-th row of the connectivity matrix to betweenness of s (follows from our formula)
	long long int numPathssz{sbd.totalSigmas[z]};
	std::queue<std::pair<VertexAppearanceSuccint, std::unordered_set<long long int>>> reverseVisit;
	// For each time that reaches the destination
	auto lastdist = 0;
	for(long long int idx_time{static_cast<long long int>(sbd.preds[z].size())-1}; idx_time >= 0; idx_time--)
	{
		//If there exists a path that reaches the vertex at that time
		if(sbd.preds[z][idx_time].second.size() > 0)
		{
			//Span each predecessor reaching the node at that time
			for (auto pred : sbd.preds[z][idx_time].second) 
			{
				auto cur_idx_time = pred.idx;
				sbd.sigmas_to_z[pred.v][cur_idx_time].second += 1;
				{
					sbd.inQueueUsages[pred.v][cur_idx_time].second += 1;
					std::unordered_set<long long int> backs{z};
					reverseVisit.push(std::make_pair(VertexAppearanceSuccint{ pred.v, pred.time, pred.idx }, backs));
					lastdist = sbd.dists[pred.v][pred.idx].second;
				}
			}
		}
	}

	bool f{true};
	// While not computing all the betweenness values
    while (!reverseVisit.empty()) {
        auto curr = reverseVisit.front();
		auto cur = curr.first;
		auto set = curr.second;
        reverseVisit.pop();
		auto thisdist = sbd.dists[cur.v][cur.idx].second;
		if(cur.v == s)
		{
		}
		else
		{
			auto cur_idx_time = cur.idx;
			auto pathstoz = sbd.sigmas_to_z[cur.v][cur_idx_time].second;
			

			if(set.count(cur.v) == 0)
			{
				currentUpdateBetweenness[cur.v] += (sbd.sigmas[cur.v][cur_idx_time].second * pathstoz / sbd.inQueueUsages[cur.v][cur_idx_time].second);
			}
			else
			{
			}
			set.emplace(cur.v);

        	for (auto pred : sbd.preds[cur.v][cur_idx_time].second) 
			{
				auto pred_idx_time = pred.idx;
				sbd.sigmas_to_z[pred.v][pred_idx_time].second += pathstoz/sbd.inQueueUsages[cur.v][cur_idx_time].second;
					
				sbd.inQueueUsages[pred.v][pred_idx_time].second += 1;
				reverseVisit.push(std::make_pair(VertexAppearanceSuccint{ pred.v, pred.time, pred.idx },set));
			}
		}
		lastdist = thisdist;
    }

	for(size_t j{0}; j < currentUpdateBetweenness.size(); j++)
	{
		if(currentUpdateBetweenness[j] > 0)
		{
			betweenness[j] += currentUpdateBetweenness[j]/static_cast<double>(numPathssz);
		}
	}

}

// Empties the stack in sbd while updating the betweenness values of the nodes of the graph
void shortestEmptyStackUpdateBetweenness(int s, ShortestBetweennessData& sbd, std::vector<double>& betweenness, std::vector<double>& fmBetweenness)
{
    // Compute the delta array values (bottom-up) and update the betweenness values
    // Subtract the sum of the s-th row of the connectivity matrix to betweenness of s (follows from our formula)
    auto connectivityCorrection = std::accumulate(sbd.totalDists.cbegin(), sbd.totalDists.cend(), int(0),
        [](auto acc, auto dist) { return acc + (dist >= 0); });
    betweenness[s] -= connectivityCorrection;
    fmBetweenness[s] -= connectivityCorrection;
    while (!sbd.stack.empty()) {
        auto cur = sbd.stack.back();
        sbd.stack.pop_back();
        // Handle "base-cases," i.e. times where the appearance is the end of a shortest (foremost) path to the vertex in question
        if (sbd.dists[cur.v][cur.time] == sbd.totalDists[cur.v])
            sbd.deltaDots[cur.v][cur.time] += static_cast<double>(sbd.sigmas[cur.v][cur.time])
                / static_cast<double>(sbd.totalSigmas[cur.v]);
        if (cur.time == sbd.foremostTimes[cur.v])
            sbd.deltaForemosts[cur.v][cur.time] += 1.0;
        for (auto pred : sbd.preds[cur.v][cur.time]) {
            // Compute the currently considered summand of the sum over successors of pred and update betweenness centrality and deltaDots
            const auto shortestSummand = (static_cast<double>(sbd.sigmas[pred.v][pred.time])
                / static_cast<double>(sbd.sigmas[cur.v][cur.time]))
                * sbd.deltaDots[cur.v][cur.time];
            sbd.deltaDots[pred.v][pred.time] += shortestSummand;
            betweenness[pred.v] += shortestSummand;
            const auto foremostSummand = (static_cast<double>(sbd.sigmas[pred.v][pred.time])
                / static_cast<double>(sbd.sigmas[cur.v][cur.time]))
                * sbd.deltaForemosts[cur.v][cur.time];
            sbd.deltaForemosts[pred.v][pred.time] += foremostSummand;
            fmBetweenness[pred.v] += foremostSummand;
        }
    }
}

void prefixDoForemostBasedSearch(const akt::Graph& g, int s, PrefixBetweennessData& pbd)
{
    auto q = std::priority_queue<akt::TemporalEdge, std::vector<akt::TemporalEdge>, decltype(&temporalEdgeGreaterTimewise)>(&temporalEdgeGreaterTimewise);
    for (auto it = g.firstEdgeFromVertex(s); (it != g.edges_cend()) && ((*it).from == s); ++it)
        q.push(*it);
    while (!q.empty()) {
        auto te = q.top();
        q.pop();
        if (pbd.foremostTimes[te.to] < 0) {
            pbd.foremostTimes[te.to] = te.when;
            pbd.stack.push_back(te.to);
            if (te.when < g.maximalTimestep())
                for (auto it = g.firstEdgeFromAppearance(te.to, te.when + 1); (it != g.edges_cend()) && ((*it).from == te.to); ++it)
                    q.push(*it);
        }
        if (te.when == pbd.foremostTimes[te.to]) {
            pbd.sigmas[te.to] += pbd.sigmas[te.from];
            pbd.preds[te.to].insert(te.from);
        }
    }
}

void prefixEmptyStackUpdateBetweenness(int s, PrefixBetweennessData& pbd, std::vector<double>& betweenness)
{
    betweenness[s] -= std::accumulate(pbd.foremostTimes.cbegin(), pbd.foremostTimes.cend(), int(0),
        [](auto acc, auto time) { return acc + (time >= 0); });
    while (!(pbd.stack.empty())) {
        auto w = pbd.stack.back();
        pbd.stack.pop_back();
        for (auto v : pbd.preds[w]) {
            // Compute the currently considered summand of the sum over successors of pred and update betweenness centrality and deltaDots
            auto summand = (static_cast<double>(pbd.sigmas[v])
                / static_cast<double>(pbd.sigmas[w]))
                * pbd.deltaDots[w];
            pbd.deltaDots[v] += summand;
            betweenness[v] += summand;
        }
    }
}

namespace akt {

    // Computes the betweenness measures
    std::pair<std::vector<double>, std::vector<double>> shortestBetweenness(const Graph& g, bool strict)
    {
		// Initialize helper structures
        auto sbd = ShortestBetweennessData(g);
        // Stores the betweenness values; initialize to 1 because of the formula for betweenness having a constant +1
        auto shortestBetweenness = std::vector<double>(g.N(), 1.0);
        auto foremostBetweenness = std::vector<double>(g.N(), 1.0);
        for (int s = 0; s < g.N(); ++s) {
            reinitializeHelperStruct(g, s, sbd, s-1);
            shortestComputeDistancesSigmas(g, strict, s, sbd);
            shortestEmptyStackUpdateBetweenness(s, sbd, shortestBetweenness, foremostBetweenness);
        }
        return { shortestBetweenness, foremostBetweenness };
    }

	// Computes the betweenness measures
    std::pair<std::vector<double>, std::vector<double>> shortestBetweennessApproxAdj(const GraphAdj& g, bool strict, int samples, int seed, double eta)
    {
		auto start1 = std::chrono::high_resolution_clock::now();
		// Initialize helper structures
        auto sbd = ShortestBetweennessDataSuccint(g);
		std::vector<int> weights(g.N(), 1);
		std::discrete_distribution<int> distribution(weights.begin(), weights.end());
		std::mt19937 generator(seed); // Seed to be set
    	std::chrono::duration<double> timeInit(0);
    	std::chrono::duration<double> timePropagate(0);
    	std::chrono::duration<double> timeBack(0);

        auto shortestBetweenness = std::vector<double>(g.N(), 0.0);
        auto foremostBetweenness = std::vector<double>(g.N(), 0.0);
		std::vector<std::vector<double>> approxesIter;
		//samples = 1;
		int prevS{-1};
		int foundpaths{0};
		for(int i=0; i < samples; i++) 
		{
        	auto shortestBetweennessInner = std::vector<double>(g.N(), 0.0);
			bool foundst = false;
			int source=-1, dest=-1;
			while(!foundst)
			{
				source = distribution(generator);
				dest = distribution(generator);
				if(source != dest)
					foundst = true;
			}
			auto start = std::chrono::high_resolution_clock::now();
            reinitializeHelperStructSuccint(g, source, sbd, prevS);
    		auto end = std::chrono::high_resolution_clock::now();
    		std::chrono::duration<double> time = end - start;
			timeInit += time;
			prevS = source;
			start = std::chrono::high_resolution_clock::now();
            bool foundpath = shortestComputeDistancesSigmasTruncatedSuccint(g, 1, source, dest, sbd);
			end = std::chrono::high_resolution_clock::now();
			time = end - start;
			timePropagate += time;
			foundpaths += foundpath? 1 : 0;
			if(foundpath)
			{
				start = std::chrono::high_resolution_clock::now();
            	updateEstimatesSuccint(source, dest, sbd, shortestBetweennessInner, g);
				end = std::chrono::high_resolution_clock::now();

				time = end - start;
				timeBack += time;
			}
				approxesIter.push_back(shortestBetweennessInner);
		}
		double max_eps = 0.0;
		double conf {eta};
		for(int i{0}; i < g.N(); i++)
		{
			std::vector<double> estim;
			for(int j{0}; j< samples; j++)
			{
				estim.push_back(approxesIter[j][i]);
			}
			double mean = std::accumulate(estim.begin(), estim.end(), 0.0, [](double l, double r) {return l+r;});
			shortestBetweenness[i] = mean;
			mean /= samples;
			double variance = 0.0;
			for(auto& es : estim)
			{
				variance += std::pow((es-mean), 2);
			}
			variance /= samples;
			double bound_eps = sqrt((2*variance*log(4*g.N()/ conf))/samples) + (7 * log(4 * g.N() / conf)/ (3*(samples-1)) );
			max_eps = std::max(max_eps, bound_eps);
		}
		double hoeffding_bound = sqrt(log(2*g.N()/conf)/(2*samples));
		std::cout << "Bound epsilon max, with " << samples << " samples is: " << max_eps << '\n';
    	std::cout << "Time to initialize structures: " << timeInit.count() <<'\n';
    	std::cout << "Time to compute forward paths: " << timePropagate.count() <<'\n';
    	std::cout << "Time to compute betweenness values: " << timeBack.count() <<'\n';
    	std::cout << "Paths to s-z found: " << foundpaths <<'\n';
		auto end1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> timeT = end1 - start1;
		std::cout << "Time needed to read and run sampling alg: " << timeT.count() << " sec\n";
        return { shortestBetweenness, foremostBetweenness };
    }

	std::pair<std::vector<double>, std::vector<double>> shortestBetweennessApproxAdjDelta(const GraphAdj& g, bool strict, int samples, int seed, int delta, double eta)
    {
		auto start1 = std::chrono::high_resolution_clock::now();
        auto sbd = ShortestBetweennessDataSuccint(g);
		std::vector<int> weights(g.N(), 1);
		std::discrete_distribution<int> distribution(weights.begin(), weights.end());
		std::mt19937 generator(seed); // Seed to be set
    	std::chrono::duration<double> timeInit(0);
    	std::chrono::duration<double> timePropagate(0);
    	std::chrono::duration<double> timeBack(0);

        auto shortestBetweenness = std::vector<double>(g.N(), 0.0);
        auto foremostBetweenness = std::vector<double>(g.N(), 0.0);
		std::vector<std::vector<double>> approxesIter;
		int prevS{-1};
		int foundpaths{0};
		for(int i=0; i < samples; i++) 
		{
        	auto shortestBetweennessInner = std::vector<double>(g.N(), 0.0);
			bool foundst = false;
			int source=-1, dest=-1;
			while(!foundst)
			{
				source = distribution(generator);
				dest = distribution(generator);
				if(source != dest)
					foundst = true;
			}
			auto start = std::chrono::high_resolution_clock::now();
            reinitializeHelperStructSuccint(g, source, sbd, prevS);
    		auto end = std::chrono::high_resolution_clock::now();
    		std::chrono::duration<double> time = end - start;
			timeInit += time;
			prevS = source;
			start = std::chrono::high_resolution_clock::now();
            bool foundpath = shortestComputeDistancesSigmasTruncatedSuccintDelta(g, 1, source, dest, sbd, delta);
			end = std::chrono::high_resolution_clock::now();
			time = end - start;
			timePropagate += time;
			foundpaths += foundpath? 1 : 0;
			if(foundpath)
			{
				start = std::chrono::high_resolution_clock::now();
            	updateEstimatesSuccintDelta(source, dest, sbd, shortestBetweennessInner, g);
				end = std::chrono::high_resolution_clock::now();

				time = end - start;
				timeBack += time;
			}
				approxesIter.push_back(shortestBetweennessInner);
		}
		double max_eps = 0.0;
		double conf {eta};
		for(int i{0}; i < g.N(); i++)
		{
			std::vector<double> estim;
			for(int j{0}; j< samples; j++)
			{
				estim.push_back(approxesIter[j][i]);
			}
			double mean = std::accumulate(estim.begin(), estim.end(), 0.0, [](double l, double r) {return l+r;});
			shortestBetweenness[i] = mean;
			mean /= samples;
			double variance = 0.0;
			for(auto& es : estim)
			{
				variance += std::pow((es-mean), 2);
			}
			variance /= samples;
			double bound_eps = sqrt((2*variance*log(4*g.N()/ conf))/samples) + (7 * log(4 * g.N() / conf)/ (3*(samples-1)) );
			max_eps = std::max(max_eps, bound_eps);
		}
		std::cout << "Bound epsilon max, with " << samples << " samples is: " << max_eps << '\n';
    	 std::cout << "Time to initialize structures: " << timeInit.count() <<'\n';
    	 std::cout << "Time to compute forward paths: " << timePropagate.count() <<'\n';
    	 std::cout << "Time to compute betweenness values: " << timeBack.count() <<'\n';
    	 std::cout << "Paths to s-z found: " << foundpaths <<'\n';
		auto end1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> timeT = end1 - start1;
		std::cout << "Time needed to read and run sampling alg: " << timeT.count() << " sec\n";
        return { shortestBetweenness, foremostBetweenness };
    }
    
} // end namespace akt
