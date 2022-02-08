#pragma once

#include "graph.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <functional>

// Possible modification: make a space-time trade-off and store the set of edges in the graph
// Possible modification: make a space-time trade-off and store the current edge of the iterator in Graph::EdgeConstIterator

namespace akt {

    // Stores data about a temporal edge
    struct TemporalEdge_64
    {
        int from;
        int to;
        long long int when;
    };

    bool temporalEdgeLessTimewise_64(const TemporalEdge_64& lhs, const TemporalEdge_64& rhs);
	// set if Edge and function pointer that acts as comparator
    using TemporalEdgeSet_64 = std::set<TemporalEdge_64, decltype(&temporalEdgeLessTimewise_64)>;

    class GraphAdj
    {
    public:
        // Helper struct for the node neighbourhood table
        // Stores all the outgoing edges of the given node at some timepoint as well as the next timestep during which the node will have outgoing edges
        struct AppearanceNeighbourhood
        {
            // Timestep of this Neighborhood, i.e., time at which the edges are out of the current node
            long long int timestep;
            // Neightbors at which the edges are going
            std::vector<int> neighbours;
        };

		/*
        // Creates a graph with (initially) no edges, containing noNodes nodes and edges whose timestamps lie in [0, maximalTimestep]
        GraphAdj(int noNodes, int maximalTimestep)
            : nodes(noNodes), edges(0), lastTime(maximalTimestep), adj(noNodes, std::vector<AppearanceNeighbourhood>(maximalTimestep + 1, AppearanceNeighbourhood{-1, std::vector<int>()}))
        { }
		*/

        // Same as Graph(int, int) except also add the edges from the set in a more efficient manner than by repeated addEdge() calls
        GraphAdj(int noNodes, long long int maximalTimestep, const TemporalEdgeSet_64& tes) : 
			//nodes(noNodes), lastTime(maximalTimestep), adj(noNodes, std::vector<AppearanceNeighbourhood>(1, AppearanceNeighbourhood{-1, std::vector<int>()})), 
			nodes(noNodes), lastTime(maximalTimestep), adj(noNodes, std::vector<AppearanceNeighbourhood>()), 
			outTimestamps(noNodes,-1)
        //    : GraphAdj(noNodes, maximalTimestep)
        {
			//this->nodes = noNodes;
			//this->lastTime = maximalTimestep;
            edges = tes.size();
			std::cout << "Tes size: " << tes.size() << '\n';
            // Add edges (without caring about the nextTimestep field for now)
			int weird{0};
            for (const auto& te : tes)
			{
				//if(adj[te.from].size() > 0)
				{
					auto pos = std::lower_bound(adj[te.from].begin(), adj[te.from].end(), te.when, [&](const AppearanceNeighbourhood& entry, long long int time){ return (entry.timestep < time) ;}) - adj[te.from].begin();
					if(pos == static_cast<int>(adj[te.from].size())) // The element is not found i.e., not in the range
						adj[te.from].push_back(AppearanceNeighbourhood{te.when, std::vector<int>{te.to}});
					else if(adj[te.from][pos].timestep > te.when) // this should never happend
						weird++;
					else
						adj[te.from][pos].neighbours.push_back(te.to);
				}
				//else
			//		adj[te.from].push_back(AppearanceNeighbourhood{te.when, std::vector<int>{te.to}});

			}
			std::cout << "Weird cases: " << weird << '\n';
			int node{0};
			for(const auto& entry : adj)
			{
				outTimestamps[node] = static_cast<long long int>(entry.size());
			/*
				std::cout << "Node: " << node << " : ";
				for(const auto& list : entry)
				{
					std::cout << "<" << list.timestep << ",(";
					for(const auto& e : list.neighbours)
						std::cout << e << ",";
					std::cout << ")>\n";
				}
			*/
				node++;

			}
        }
		
        // Returns the number of nodes in the graph
        int N() const { return nodes; }

        // Returns the number of temporal edges in the graph
        long long int M() const { return edges; }

        // Returns the number of timesteps in the graph (so maximalTimesteps() + 1)
        long long int T() const { return lastTime + 1; }

        // Returns the last timestep in the graph
        long long int maximalTimestep() const { return lastTime; }

        // Returns the O(nT) adjacency list for the graph
        const std::vector<std::vector<AppearanceNeighbourhood>>& adjacencyList() const { return std::cref(adj); }

        // returns a vector contaiing the number of out timesteps for each vertex
        const std::vector<long long int>& getOutTimestamps() const { return std::cref(outTimestamps); }

		/*
        // Helper class which lets the user iterate over the edges
        class EdgeConstIterator
        {
			// Graph can access private members
            friend class GraphAdj;
        public:
            EdgeConstIterator& operator++()
            {
                // Check for simple cases
                if ((nodeNo == -1) || (timeNo == -1))
                    return *this;
                if ((static_cast<uint64_t>(edgeNo) + 1) < ((*data)[nodeNo][timeNo].neighbours.size())) {
                    ++edgeNo;
                    return *this;
                }
                // Increment the iterator for the complicated cases
                goToFollowingNonemptyCell();
                return *this;
            }

            TemporalEdge operator*() const
            {
                return TemporalEdge{ nodeNo, (*data)[nodeNo][timeNo].neighbours[edgeNo], timeNo };
            }

            bool operator==(const EdgeConstIterator& rhs)
            {
                return (this->data == rhs.data)
                    && (this->nodeNo == rhs.nodeNo)
                    && (this->timeNo == rhs.timeNo)
                    && (this->edgeNo == rhs.edgeNo);
            }

            bool operator!=(const EdgeConstIterator& rhs)
            {
                return !(*this == rhs);
            }


        private:
            EdgeConstIterator(const std::vector<std::vector<AppearanceNeighbourhood>>* data, int nodeNo, int timeNo, int edgeNo)
                :data(data), nodeNo(nodeNo), timeNo(timeNo), edgeNo(edgeNo) { }

            // Finds the first non-empty cell following the current one (or sets all the stuff to -1 if there is no such cell)
            void goToFollowingNonemptyCell()
            {
                if ((nodeNo < 0) || (timeNo < 0))
                    return;
                // See if the current node has any further temporal edges
                if ((*data)[nodeNo][timeNo].nextTimestep >= 0) {
                    timeNo = (*data)[nodeNo][timeNo].nextTimestep;
                    edgeNo = 0;
                    return;
                }
                // Find the first node with any temporal edges
                ++nodeNo;
                while ((nodeNo < ((*data).size())) && ((*data)[nodeNo][0].neighbours.empty()) && ((*data)[nodeNo][0].nextTimestep < 0))
                    ++nodeNo;
                if (nodeNo >= (*data).size())
                    nodeNo = timeNo = edgeNo = -1;
                else if ((*data)[nodeNo][0].neighbours.empty()) {
                    timeNo = (*data)[nodeNo][0].nextTimestep;
                    edgeNo = 0;
                }
                else {
                    timeNo = edgeNo = 0;
                }
            }

            // The data where the set of edges is taken from
            const std::vector<std::vector<AppearanceNeighbourhood>>* const data;
            // The position within the edge set
            int nodeNo, timeNo, edgeNo;
        };

        // Returns an iterator to the beginning of the edge set
        EdgeConstIterator edges_cbegin() const
        {
            if (nodes <= 0)
                return EdgeConstIterator(&adj, -1, -1, -1);
            auto res = EdgeConstIterator(&adj, 0, 0, 0);
            if (adj[0][0].neighbours.size() == 0)
                res.goToFollowingNonemptyCell();
            return res;
        }

        // Returns an iterator to the first edge starting at the given vertex (with time 0)
        EdgeConstIterator firstEdgeFromVertex(int v) const
        {
            return firstEdgeFromAppearance(v, 0);
        }

        // Returns an iterator to the first edge starting at the given vertex appearance
        EdgeConstIterator firstEdgeFromAppearance(int v, int time) const
        {
            auto res = EdgeConstIterator(&adj, v, time, 0);
            if (adj[v][time].neighbours.size() == 0)
                res.goToFollowingNonemptyCell();
            return res;
        }

        EdgeConstIterator edges_cend() const
        {
            return EdgeConstIterator(&adj, -1, -1, -1);
        }
		*/

    private:
        // The number of nodes in the graph
        int nodes;
        // The total number of temporal edges in the graph
        long long int edges;
        // The "timespan" of the graph, i.e. the timestamps on the temporal edges are in [0, lastTime]
        long long int lastTime;
        // The huge node lookup table (O(lastTime) for each node in the graph, so O(nT) total)
        std::vector<std::vector<AppearanceNeighbourhood>> adj;
        // for each vertex keeps the number of different entries in adj, i.e., number of different timesteps 
		// from which there exists at least one outgoing edge from the vertex + 1 (to init structures)
		std::vector<long long int> outTimestamps;
    };

    // Reads a graph from the given input stream; automatically performs a basic reduction (removes times with no temporal edges within them)
    // The proper file format is a sequence of lines 
    //  with each line starting with two (integer or string) nodeIDs and the (integer) timepoint at which the edge appears
    // Self-loops, duplicate edges, or invalid lines are ignored
    // Simple examples of a correct graph would be
    // 1 2 0
    // 1 3 0
    // 2 4 1
    // 3 4 2
    // or
    // a 13 3
    // a 13 1
    // 13 f5 2
    // a f5 0
    // If directed == false, then the symmetric edge is added for every edge
    // Returns the Graph read alongside a table of mappings between the assigned node IDs to the original node IDs from the input
    //std::pair<Graph, std::vector<std::string>> readReduceGraph(std::istream& is, bool directed = false);
    std::pair<GraphAdj, std::vector<std::string>> readReduceGraphToAdj(std::istream& is, bool directed=false);

} // end namespace akt
