# ONBRA Rigorous Estimation of the Temporal Betweenness Centrality in Temporal Temporal Networks (Santoro and Sarpe, TheWebConf 2022) 

## How to build

- `mkdir build.build`
- `cd build.build`
- `cmake ../src`
- `make`

## How to use

 - Build (see `How to build` above)
 - Use the `./onbra` executable by supplying a graph in the proper format (described below)
 - Use `./onbra -h` to see a list of available flags and options
 - Note: option `-E` selects if you need to compute the Temporal Betweenness for shortest paths (set it to `1`) or for shortest $\delta$-restless walks (set it to `3` and in such case you will also nead the parameter `-D`)
  
 Example to use ONBRA for shortest paths with a sample size of 1000 nodes over 10 executions and appending all the results in "result.txt":
 - ./onbra -f <filename> -d -s -E 1 -S 1000 -I 10 &> result.txt
 Example to use ONBRA for restless walks with a sample size of 1000 nodes, delta 3200, over 10 executions and appending all the results in "result.txt":
 - ./onbra -f <filename> -d -s -E 3 -S 1000 -I 10 -D 3200 &> result.txt

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
