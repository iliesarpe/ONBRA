#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <tuple>
#include <future>

#include <boost/program_options.hpp>

#include "adjList.h"
#include "graph.h"
#include "algorithms.h"
#include "ctpl.h"

namespace po = boost::program_options;

struct BenchmarkResults 
{
    std::vector<double> shortest, foremost, strictShortest, strictForemost, prefix;
    int n;
    std::vector<std::string> inputIds;
    double nonStrictTime = -1.0;
    double strictTime = -1.0;
    double prefixTime = -1.0;
};

struct BenchmarkSettings
{
    bool runStrict = true;
    bool runNonStrict = true;
    bool runPrefix = true;
    bool edgesDirected = false;
    bool originalNodeIds = false;
    bool readFromFile = false;
    std::string filename;
};



auto readGraphFromFile(const BenchmarkSettings& bs)
{
    std::ifstream ifs{ bs.filename };
    if (!ifs) {
        std::cout << "Error trying to open the file \"" << bs.filename << "\"\n";
        throw 1;
    }
    return akt::readReduceGraph(ifs, bs.edgesDirected);
}

BenchmarkResults readGraphRunShortestBetweenness(const BenchmarkSettings& bs, int samples, int eoa, int execs=1, int delta=0, double eta=0.1)
{
	if(eoa==1)
	{
		// Not the perfect solution from the perspective of error handling (opening the file and not checking for success), but simple to write cleanly)
		auto start = std::chrono::high_resolution_clock::now();
    	std::ifstream ifs{ bs.filename };
    	auto [g, ids] = akt::readReduceGraphToAdj(ifs, bs.edgesDirected);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time = end - start;
		std::cout << "Time needed to read and preprocess data: " << time.count() << " sec\n";

		std::clog << "Graph read: " << g.N() << " nodes, ";
		if (bs.edgesDirected)
			std::clog << g.M() << " (unique) directed edges, ";
		else
			std::clog << g.M() / 2 << " (unique) edges, ";
		std::clog << g.T() << " (non-empty) timesteps\n";

		BenchmarkResults res;

		int seeds[10] {10, 254, 4321, 91283, 39483, 10239102, 11864, 65683, 891273, 7684};
		std::vector<std::vector<double>> approxes(g.N());
		for(int j{0}; j<execs; j++)
		{
			start = std::chrono::high_resolution_clock::now();
			std::cout << "Samples used: " << samples << '\n';
			std::tie(res.strictShortest, res.strictForemost) = shortestBetweennessApproxAdj(g, true, samples, seeds[j], eta);
			end = std::chrono::high_resolution_clock::now();
			time = end - start;
			std::cout << "Time needed to read and run sampling alg: " << time.count() << " sec\n";

			res.inputIds = std::move(ids);
			res.n = g.N();
			for (int i = 0; i < res.n; ++i) {
					//std::cout << "Node " << i << ": " << res.strictShortest[i]/(1.*samples) << "\n";
					approxes[i].push_back(res.strictShortest[i]/(1.*samples));
			}
		}
		for(size_t j{0}; j<approxes.size(); j++)
		{
			std::cout << "Node " << j << ": ";
			for(auto& approx : approxes[j])
			{
				std::cout << approx << " ";
			}
			std::cout << '\n';
		}	
    	return res;
	}
	else if(eoa==2) // Parallel Iterations - Sequential Sampling Succint
	{
		// Not the perfect solution from the perspective of error handling (opening the file and not checking for success), but simple to write cleanly)
		auto start = std::chrono::high_resolution_clock::now();
    	std::ifstream ifs{ bs.filename };

    	auto [g, ids] = akt::readReduceGraphToAdj(ifs, bs.edgesDirected);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time = end - start;
		std::cout << "Time needed to read and preprocess data: " << time.count() << " sec\n";

		std::clog << "Graph read: " << g.N() << " nodes, ";
		if (bs.edgesDirected)
			std::clog << g.M() << " (unique) directed edges, ";
		else
			std::clog << g.M() / 2 << " (unique) edges, ";
		std::clog << g.T() << " (non-empty) timesteps\n";

		BenchmarkResults res;

		int seeds[10] {10, 254, 4321, 91283, 39483, 10239102, 11864, 65683, 891273, 7684};
		std::vector<std::vector<double>> approxes(g.N());
		//Define lambda function
		auto exectit{[&](int ID, int samp, int seed){return shortestBetweennessApproxAdj(g, true, samp, seed, eta);} };
		std::vector<std::future<std::pair<std::vector<double>, std::vector<double>>>> futures;
		ctpl::thread_pool p(execs);

		for(int j{0}; j<execs; j++)
		{
			start = std::chrono::high_resolution_clock::now();
			std::cout << "Samples used: " << samples << '\n';
			//std::tie(res.strictShortest, res.strictForemost) = shortestBetweennessApproxAdj(g, true, samples, seeds[j]);
			futures.push_back(p.push(std::ref(exectit), samples, seeds[j]));
			end = std::chrono::high_resolution_clock::now();
			time = end - start;
			//std::cout << "Time needed to read and run sampling alg: " << time.count() << " sec\n";

			/*
			res.inputIds = std::move(ids);
			res.n = g.N();
			for (int i = 0; i < res.n; ++i) {
					//std::cout << "Node " << i << ": " << res.strictShortest[i]/(1.*samples) << "\n";
					approxes[i].push_back(res.strictShortest[i]/(1.*samples));
			}
			*/
		}
		for(int j{0}; j<execs; j++)
		{
			std::tie(res.strictShortest, res.strictForemost) = futures[j].get();
			res.inputIds = std::move(ids);
			res.n = g.N();
			for (int i = 0; i < res.n; ++i) {
					//std::cout << "Node " << i << ": " << res.strictShortest[i]/(1.*samples) << "\n";
					approxes[i].push_back(res.strictShortest[i]/(1.*samples));
			}
		}
		for(size_t j{0}; j<approxes.size(); j++)
		{
			std::cout << "Node " << j << ": ";
			for(auto& approx : approxes[j])
			{
				std::cout << approx << " ";
			}
			std::cout << '\n';
		}	
    	return res;
	}
	else if(eoa==3)
	{
		// Not the perfect solution from the perspective of error handling (opening the file and not checking for success), but simple to write cleanly)
		auto start = std::chrono::high_resolution_clock::now();
    	std::ifstream ifs{ bs.filename };
    	auto [g, ids] = akt::readReduceGraphToAdj(ifs, bs.edgesDirected);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time = end - start;
		std::cout << "Time needed to read and preprocess data: " << time.count() << " sec\n";

		std::clog << "Graph read: " << g.N() << " nodes, ";
		if (bs.edgesDirected)
			std::clog << g.M() << " (unique) directed edges, ";
		else
			std::clog << g.M() / 2 << " (unique) edges, ";
		std::clog << g.T() << " (non-empty) timesteps\n";

		BenchmarkResults res;

		int seeds[10] {10, 254, 4321, 91283, 39483, 10239102, 11864, 65683, 891273, 7684};
		std::vector<std::vector<double>> approxes(g.N());
		for(int j{0}; j<execs; j++)
		{
			start = std::chrono::high_resolution_clock::now();
			std::cout << "Samples used: " << samples << '\n';
			std::tie(res.strictShortest, res.strictForemost) = shortestBetweennessApproxAdjDelta(g, true, samples, seeds[j], delta, eta);
			end = std::chrono::high_resolution_clock::now();
			time = end - start;
			std::cout << "Time needed to read and run sampling alg: " << time.count() << " sec\n";

			res.inputIds = std::move(ids);
			res.n = g.N();
			for (int i = 0; i < res.n; ++i) {
					//std::cout << "Node " << i << ": " << res.strictShortest[i]/(1.*samples) << "\n";
					approxes[i].push_back(res.strictShortest[i]/(1.*samples));
			}
		}
		for(size_t j{0}; j<approxes.size(); j++)
		{
			std::cout << "Node " << j << ": ";
			for(auto& approx : approxes[j])
			{
				std::cout << approx << " ";
			}
			std::cout << '\n';
		}	
    	return res;
	}
}

void outputBenchmarkResults(const BenchmarkSettings& bs, const BenchmarkResults& br)
{
    std::clog << "Time for algorithms (in seconds):\n";
    std::clog << "Non-strict, strict, prefix foremost\n";
    std::clog << br.nonStrictTime << ", " << br.strictTime << ", " << br.prefixTime << '\n';
    std::cout << br.nonStrictTime << ", " << br.strictTime << ", " << br.prefixTime << '\n';
    std::clog << "Computed betweenness measures:\n";
    std::clog << "Node, non-strict shortest, non-strict shortest foremost, strict shortest, strict shortest foremost, prefix betweenness\n";
    for (int i = 0; i < br.n; ++i) {
        if (bs.originalNodeIds)
            std::cout << br.inputIds[i];
        else
            std::cout << i;
        std::cout << ", ";
        if (bs.runNonStrict)
            std::cout << br.shortest[i] << ", " << br.foremost[i] << ", ";
        else
            std::cout << "-1, -1, ";
        if (bs.runStrict)
            std::cout << br.strictShortest[i] << ", " << br.strictForemost[i] << ", ";
        else
            std::cout << "-1, -1, ";
        if (bs.runPrefix)
            std::cout << br.prefix[i] << '\n';
        else
            std::cout << "-1\n";
    }
}

int main (int argc, char** argv)
{
    BenchmarkSettings bs;
    po::options_description desc("usage: btwBenchmark [options]\nRuns the different betweenness centrality algorithms on the graph input via stdin (or file if -f used). Available options");
    desc.add_options()
		("help,h", "write help message")
		("filename,f", po::value<std::string>(&(bs.filename)), "instead of reading the graph from stdin, use the file given in the argument")
        ("graph-directed,d", "interpret the edges in the graph as directed edges")
		//("no-strict,s", "don't run the strict shortest (foremost) betweenness algorithm")
		//("no-non-strict,n", "don't run the non-strict shortest (foremost) betweenness algorithm")
		//("no-prefix,p", "don't run the the strict prefix-foremost betweenness algorithm")
		//("original-node-ids,u", "in the output, use original node ids instead of the arbitrary integer ids assigned by the program")
		("exact-or-approximate,E", po::value<int>(), "Select the mode:\n 1) Execute ONBRA on STPs with all the <I> iterations sequentially\n 2) Execute ONBRA on STPs with all the <I> iterations in parallel (each iteration is still sequential!)\n 3) Execute ONBRA on RTPs with all the <I> iterations in parallel (each iteration is still sequential!)")
		("iterations,I", po::value<int>(), "Iterations to execute ONBRA  (with different seeds)")
		("samples,S",po::value<int>(), "Number of vertex pairs to be selected")
		("prob,T",po::value<int>(), "Probability eta")
		("delta,D",po::value<int>(), "Bound on restless duration");
        
    po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << '\n';
        return 0;
    }
    bs.readFromFile = vm.count("filename") > 0;
    bs.edgesDirected = vm.count("graph-directed") > 0;
    bs.runStrict = vm.count("no-strict") <= 0;
    bs.runNonStrict = vm.count("no-non-strict") <= 0;
    bs.runPrefix = vm.count("no-prefix") <= 0;
    bs.originalNodeIds = vm.count("originalNodeIds") > 0;
	int s { 0};
	double eta{ 0.1 };
	int delt{ 0};
	int execs{1};
	if(vm.count("samples"))
	{
		s += vm["samples"].as<int>();
	}
	if(vm.count("delta"))
	{
		delt += vm["delta"].as<int>();
	}
	int eoa{0};
	if(vm.count("exact-or-approximate"))
	{
		eoa += vm["exact-or-approximate"].as<int>();
	}
	if(vm.count("iterations"))
	{
		execs = vm["iterations"].as<int>();
	}
	if(vm.count("prob"))
	{
		s += vm["prob"].as<double>();
	}
    
    try {
        auto br = readGraphRunShortestBetweenness(bs, s, eoa, execs, delt, eta);
    } catch (...) {
        return 1;
    }
}
