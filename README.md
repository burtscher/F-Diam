# F-Diam

## Description

F-Diam is a C++/OpenMP code for quickly computing the exact diameter of large sparse graphs. It utilizes techniques such as Winnowing to minimize the number of needed BFS calls. Each key level-synchronous BFS is run in parallel.

## Installation

Follow these instructions to install, compile, and run the code.

1. Clone the code using "git clone https://github.com/burtscher/F-Diam.git"
2. Navigate to the F-Diam directory with "cd F-Diam"
3. Execute "make"
4. Run "./fdiam graphs/internet.egr"
5. Optionally, download more inputs using "bash setup_full.sh"

To convert your own graphs into our format, follow the instructions at https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/.
