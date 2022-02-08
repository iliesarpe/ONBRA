#pragma once

#include <iostream>
#include <set>
#include <string>
#include <vector>

// Possible modification: make a space-time trade-off and store the set of edges in the graph
// Possible modification: make a space-time trade-off and store the current edge of the iterator in Graph::EdgeConstIterator

namespace akt {

    // Stores data about a temporal edge
    struct TemporalEdge
    {
        int from;
        int to;
        int when;
    };

    bool temporalEdgeLessTimewise(const TemporalEdge& lhs, const TemporalEdge& rhs);
	// set if Edge and function pointer that acts as comparator
    using TemporalEdgeSet = std::set<TemporalEdge, decltype(&temporalEdgeLessTimewise)>;

    class Graph
    {
    public:
        // Helper struct for the node neighbourhood table
        // Stores all the outgoing edges of the given node at some timepoint as well as the next timestep during which the node will have outgoing edges
        struct AppearanceNeighbourhood
        {
            // Outgoing edges from the node
            std::vector<int> neighbours;
            // Next timestep where the node will have outgoing edges (or -1 if no such timestep exists)
            int nextTimestep;
        };

        // Creates a graph with (initially) no edges, containing noNodes nodes and edges whose timestamps lie in [0, maximalTimestep]
        Graph(int noNodes, int maximalTimestep)
            : nodes(noNodes), edges(0), lastTime(maximalTimestep), adj(noNodes, std::vector<AppearanceNeighbourhood>(maximalTimestep + 1, AppearanceNeighbourhood{ std::vector<int>(), -1 }))
        { }

        // Same as Graph(int, int) except also add the edges from the set in a more efficient manner than by repeated addEdge() calls
        Graph(int noNodes, int maximalTimestep, const TemporalEdgeSet& tes)
            : Graph(noNodes, maximalTimestep)
        {
            edges = tes.size();
            // Add edges (without caring about the nextTimestep field for now)
            for (const auto& te : tes)
                adj[te.from][te.when].neighbours.push_back(te.to);
            // Compute the nextTimestep values for all nodes at all times
            for (int i = 0; i < noNodes; ++i) {
                // Initialize lastNonEmpty to the last timestep where node i has some outgoing edges
                int lastNonEmpty = maximalTimestep;
                for (; (lastNonEmpty >= 0) && (adj[i][lastNonEmpty].neighbours.size() == 0); --lastNonEmpty)
                    ;
                for (int curTime = lastNonEmpty - 1; curTime >= 0; ) {
                    for (; (curTime >= 0) && (adj[i][curTime].neighbours.size() == 0); --curTime)
                        adj[i][curTime].nextTimestep = lastNonEmpty;
                    if (curTime >= 0) {
                        adj[i][curTime].nextTimestep = lastNonEmpty;
                        lastNonEmpty = curTime;
                        --curTime;
                    }
                }
            }
        }

        // Adds the temporal edge to the graph
        // NOTE: following conditions must be satisfied for successful addition:
        //  - Both source and target nodes must exist in the graph
        //  - The timestamp on the edge must be less than or equal to current maximal timestamp--otherwise you need to extend the current max timestamp first
        //      (The extension isn't done automatically because of its high cost) 
        // Returns true on successful addition of the edge and false otherwise
        // Note: addEdge() doesn't check for duplicates so the same edge can be added to the graph multiple times
        bool addEdge(const TemporalEdge& te)
        {
            if ((te.from >= nodes) || (te.to >= nodes) || (te.when > lastTime))
                return false;

            // Add edge
            adj[te.from][te.when].neighbours.push_back(te.to);
            ++edges;
            // Update next non-empty cell, if needed
            if (adj[te.from][te.when].neighbours.size() == 1) {
                int time = te.when - 1;
                for (; (time >= 0) && (adj[te.from][time].neighbours.size() == 0); --time)
                    adj[te.from][time].nextTimestep = te.when;
                if (time >= 0)
                    adj[te.from][time].nextTimestep = te.when;
            }

            return true;
        }

        // Returns the number of nodes in the graph
        int N() const { return nodes; }

        // Returns the number of temporal edges in the graph
        int M() const { return edges; }

        // Returns the number of timesteps in the graph (so maximalTimesteps() + 1)
        int T() const { return lastTime + 1; }

        // Returns the last timestep in the graph
        int maximalTimestep() const { return lastTime; }

        // Returns the O(nT) adjacency list for the graph
        const std::vector<std::vector<AppearanceNeighbourhood>>& adjacencyList() const { return adj; }

        // Helper class which lets the user iterate over the edges
        class EdgeConstIterator
        {
			// Graph can access private members
            friend class Graph;
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
                while ((nodeNo < static_cast<int>(((*data).size()))) && ((*data)[nodeNo][0].neighbours.empty()) && ((*data)[nodeNo][0].nextTimestep < 0))
                    ++nodeNo;
                if (nodeNo >= static_cast<int>((*data).size()))
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

    private:
        // The number of nodes in the graph
        int nodes;
        // The total number of temporal edges in the graph
        int edges;
        // The "timespan" of the graph, i.e. the timestamps on the temporal edges are in [0, lastTime]
        int lastTime;
        // The huge node lookup table (O(lastTime) for each node in the graph, so O(nT) total)
        std::vector<std::vector<AppearanceNeighbourhood>> adj;
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
    std::pair<Graph, std::vector<std::string>> readReduceGraph(std::istream& is, bool directed = false);

} // end namespace akt
