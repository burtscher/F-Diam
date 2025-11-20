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

## Publication
If you use F-Diam in your work, please cite the following publication: <br>
Cameron Bradley, Anju Mongandampulath Akathoott, and Martin Burtscher. Fast Exact Diameter Computation of Sparse Graphs. Proceedings of the 54th International Conference on Parallel Processing. September 2025. 
[[doi]](https://doi.org/10.1145/3754598.3754620) [[paper]](https://userweb.cs.txstate.edu/~burtscher/papers/icpp25.pdf) [[slides]](https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fuserweb.cs.txstate.edu%2F~burtscher%2Fslides%2Ficpp25.pptx&wdOrigin=BROWSELINK)

*This work has been supported by the National Science Foundation under Award Number 1955367.*
