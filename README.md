# ONBRA Rigorous Estimation of the Temporal Betweenness Centrality in Temporal Networks (Santoro and Sarpe, TheWebConf 2022) 

## How to build

- `mkdir build.build`
- `cd build.build`
- `cmake ../src`
- `make`

## How to use

 - Build (see `How to build` above)
 - Use the `./onbra` executable by supplying a graph in the proper format (described below)
 - Use `./onbra -h` to see a list of available flags and options
 - Note: option `-E` selects if you need to compute the Temporal Betweenness for shortest paths (set it to `1`) or for shortest $\delta$-restless walks (set it to `3` and in such case you will also need the parameter `-D`)
  
 Example to use ONBRA for shortest paths with a sample size of 1000 pairs of nodes over 10 executions and appending all the results in "result.txt":
`./onbra -f <filename> -d -s -E 1 -S 1000 -I 10 &> result.txt`
 
 Example to use ONBRA for restless walks with a sample size of 1000 pairs of nodes, delta 3200, over 10 executions and appending all the results in "result.txt":
`./onbra -f <filename> -d -s -E 3 -S 1000 -I 10 -D 3200 &> result.txt`

## Output of ONBRA
A sample of the output that you may obtain by running ONBRA is the following (we comment each line by adding an arrow (&rarr;)):

Samples used: $S$ &rarr; $S$: is the sample size you provided in input. \
Bound epsilon, max with $S$ samples is: $\varepsilon$ &rarr; $\varepsilon$: is a bound on the supremum deviation in current iteration using the empirical Bernstein bound (see paper). \
Time to initizialize structures: $t_1$  &rarr; $t_1$: time needed to initialize internal structures to ONBRA \
Time to compute forward paths: $t_2$ &rarr; $t_2$: time to compute paths for sampled pairs of nodes \
Time to compute betweenness values: $t_3$ &rarr; $t_3$: time used to process identified paths to update node values \
Paths to s-z found: $P$ &rarr; $P$: number of pairs of nodes $(s,z)$ sampled for which there exists at least a path from $s$ to $z$. \
Time needed to read and run sampling alg: $t_{tot}$ &rarr; $t_{tot}$: total time to run an iteration of ONBRA (is at least $t_1+t_2+t_3$) \
Node 0: $[b(0)_1' \cdots b(0)_I']$ -> each $b(0)_i'$ is the estimate obtained by ONBRA in $i$-th iteration ($i \in [1,I]$) for node 0 \
$\cdots$ \
Node $n-1$: $[b(n-1)_1' \cdots b(n-1)_I']$ -> each $b(n-1)_i'$ is the estimate obtained by ONBRA in $i$-th iteration ($i \in [1,I]$) for node $n-1$

## Graph format for input into the ONBRA

Temporal graphs which are read by the benchmark suite need to have the following form:
 - A graph is represented by a sequence of lines, with each corresponding to a temporal edge in the graph
 - Node IDs should be non negative integers (starting from 0 included), such as `42` or `302` are valid node IDs. 
 - Timestamps must be positive integers (strictly greater than 0) but such that they fit inside 64-bit signed integer
 - Each line of the input must start with the following description of an edge: ID of the origin (tail) node, ID of the destination (head) node, timestamp, all separated by (non-newline) whitespace.
 - We assume the input network is pre-processed such that edges appear time-ordered, and nodes appear sequentially, i.e., node id $i$ cannot appear on one edge before $i-1$ has not been seen on some other edge. We provide a script to preprocess a dataset in the folder `utils`.
 - Self-loops and duplicate edges are allowed, however they will be ignored
Example of a valid temporal network:\
    0 1 1\
    1 2 1\
    1 2 2\
    2 3 3
