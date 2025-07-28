#!/bin/bash

mkdir -p graphs
cd graphs

wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/2d-2e20.sym.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/amazon0601.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/as-skitter.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/citationCiteseer.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/cit-Patents.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/coPapersDBLP.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/delaunay_n24.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/europe_osm.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/in-2004.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/kron_g500-logn21.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/rmat16.sym.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/rmat22.sym.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/soc-LiveJournal1.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/USA-road-d.NY.egr
wget https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/USA-road-d.USA.egr

cd ..

echo "Done"
