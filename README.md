# Dynamic programming for the time-dependent traveling salesman problem with time windows
Source code to replicate the experiments from the article.

## Abstract
The Time-Dependent Traveling Salesman Problem with Time Windows (TDTSPTW) is a variant of the well-known Traveling Salesman Problem with Time Windows (TSPTW) in which travel times are not assumed to be constant. The TDTSPTW accounts the effects of congestion at the planning level, being particularly suited for distribution problems in large cities. In this article we develop a labeling-based algorithm for the TDTSPTW that incorporates state-of-the-art components from the related literature.  We propose a new state-space relaxation specifically designed for the time-dependent context. Extensive computational experiments show the effectiveness of the overall approach and the impact of the new relaxation, outpeforming several recent algorithms proposed for the TDTSPTW. In addition, we provide evidence showing that our approach also improves the recent results reported for the Minimum Tour Duration Problem.

## Getting started
The following instructions will guide you through the steps to execute the experiments from the article.

### Prerequisites
- Python >= 3.6 [(more info)](https://www.python.org/)
- CPLEX >= 12.8 [(more info)](https://www.ibm.com/products/ilog-cplex-optimization-studio)
- Boost Graph Library >=1.66 [(more info)](https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)
    - On Linux: ```sudo apt-get install libboost-all-dev```
- CMake >= 2.8.4 [(more info)](https://cmake.org/)
    - On Linux: ```sudo apt-get install cmake```
- C++14 or higher [(more info)](https://es.wikipedia.org/wiki/C%2B%2B14)

### Built with
- Kaleidoscope: A tool to visualize the outputs of Optimization Problems [(more info)](https://github.com/gleraromero/kaleidoscope)
- Runner: A script to ease the process of running experiments [(more info)](https://github.com/gleraromero/runner)
- GOC lib: A library that includes interfaces for using (Mixed Integer) Linear Programming solvers, and some useful resources [(more info)](https://github.com/gleraromero/goc).

### Running the experiments.
1. Add environment variables with the paths to the libraries.
    1. Add two environment variables to bash with CPLEX include and library paths.
        1. ```export CPLEX_INCLUDE=<path_to_cplex_include_dir>```
            - Usually on Linux: _/opt/ibm/ILOG/CPLEX_Studio\<VERSION\>/cplex/include_
        1. ```export CPLEX_BIN=<path_to_cplex_lib_binary_file>```
            - Usually on Linux: _/opt/ibm/ILOG/CPLEX_Studio\<VERSION\>/cplex/lib/x86-64_linux/static_pic/libcplex.a_
    1. Add two environment variables to bash with BOOST Graph Library include and library paths.
        1. ```export BOOST_INCLUDE=<path_to_boost_include_dir>```
            - Usually on Linux: _/usr/include_
        1. ```export BOOST_BIN=<path_to_boost_lib_binary_file>```
            - Usually on Linux: _/usr/lib/x86_64-linux-gnu/libboost_graph.a_
1. Go to the ejor2019 root directory.
1. Execute ```python3 runner/runner.py <experiment_file>```
1. The execution output will be continually saved to the output folder.

> Experiment files are located in the _experiments_ folder. For more information see Section [Experiments](#Experiments)

### Experiments
There are eight experiment files. Which together comprise all the experiments carried out in the article.
* _Section 6.1_: arigliano_reduced.json guerriero.json
* _Section 6.2_: arigliano_reduced.json arigliano.json
* _Section 6.3_: adamo.json
* _Section 6.4_: ascheuer.json gendreau.json potvin-bengio.json olhmann-thomas.json
* _Section 6.5_: tdcespp-root.json tdop-root.json tdptptwpd-root.json tdtsptw-root.json

### Visualizing the experiment results.
1. Go to folder Kaleidoscope and run the index.html on a web browser.
1. Add the output file.
1. Select the experiments.
1. Add some attributes to visualize.
1. Click on Refresh.
1. If more details on an experiment are desired click on the + icon in a specific row.

## Built With
* [JSON for Modern C++](https://github.com/nlohmann/json)
* [Boost Graph Library](https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)

## Authors
- Gonzalo Lera-Romero
- Juan José Miranda-Bront
- Francisco Soulignac

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
