/*
F-Diam is a C++/OpenMP code for quickly computing the exact diameter of large sparse graphs.

Copyright (c) 2025, Cameron Bradley, Anju Mongandampulath Akathoott, and Martin Burtscher

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

URL: The latest version of this code is available at https://github.com/burtscher/F-Diam.git.

Publication: This work is described in detail in the following paper.
Cameron Bradley, Anju Mongandampulath Akathoott, and Martin Burtscher. "Fast Exact Diameter Computation of Sparse Graphs." Proceedings of the 54th International Conference on Parallel Processing. September 2025.

Sponsor: This work has been supported by the National Science Foundation under Award #1955367.
*/


#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <climits>
#include <omp.h>
#include <sys/time.h>
#include <vector>
#include "ECLgraph.h"


static ECLgraph g;
static int iter = 0;
static int* label;
static int* dist;
static int* wl1;
static int* wl2;
static int* wl_index;


struct CPUTimer
{
  timeval beg, end;
  CPUTimer() {}
  ~CPUTimer() {}
  void start() {gettimeofday(&beg, NULL);}
  double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};


static int distance(const int v, int& sum, int diameter)
{
  // level-by-level BFS
  iter++;
  label[v] = iter;
  int g_level = 0;
  int size1 = g.nindex[v + 1] - g.nindex[v];
  sum = size1 + 1;
  wl_index[0] = 0;

  #pragma omp parallel default(none) shared(wl1, g, label, iter, size1, wl_index, v) reduction(max: g_level) reduction(+: sum)
  {
    int level = 0;
    int l_sum = 0;
    int curVal;
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    std::vector<int> wl2_v;
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
    float ratio;
    level++;

    // for loop through first adjacency list
    #pragma omp for
    for (int u = beg; u < end; u++) {
      const int n = g.nlist[u];
      label[n] = iter;
      wl1[u - beg] = n;
    }

    do {
      level++;
      ratio = static_cast<float>(size1) / g.nodes;
      if (ratio < 0.1) {
        long start = tid * static_cast<long>(size1) / num_threads;
        long end = (tid + 1) * static_cast<long>(size1) / num_threads;
        for (int i = start; i < end; i++) {
          const int u = wl1[i];
          for (int j = g.nindex[u]; j < g.nindex[u + 1]; j++) {
            const int n = g.nlist[j];
            #pragma omp atomic capture
            {
              curVal = label[n];
              label[n] = iter;
            }
            if (curVal != iter) {
              wl2_v.push_back(n);
            }
          }
        }
        l_sum += wl2_v.size();
        wl_index[tid + 1] = wl2_v.size();
        #pragma omp barrier
        #pragma omp single
        {
          for (int i = 2; i <= num_threads; i++) {
            wl_index[i] += wl_index[i - 1];
          }
          size1 = wl_index[num_threads];
        }
        int counter = 0;
        for (int i = wl_index[tid]; i < wl_index[tid + 1]; i++) {
          wl1[i] = wl2_v[counter];
          counter++;
        }
      } else {
        long start = tid * static_cast<long>(g.nodes) / num_threads;
        long end = (tid + 1) * static_cast<long>(g.nodes) / num_threads;
        for (int i = start; i < end; i++) {
          if (label[i] != iter) {
            for (int j = g.nindex[i]; j < g.nindex[i + 1]; j++) {
              const int n = g.nlist[j];
              if (label[n] == iter) {
                wl2_v.push_back(i);
                break;
              }
            }
          }
        }
        l_sum += wl2_v.size();
        wl_index[tid + 1] = wl2_v.size();
        #pragma omp barrier
        #pragma omp single
        {
          {
            for (int i = 2; i <= num_threads; i++) {
              wl_index[i] += wl_index[i - 1];
            }
            size1 = wl_index[num_threads];
          }
        }
        int counter = 0;
        for (int i = wl_index[tid]; i < wl_index[tid + 1]; i++) {
          const int u = wl2_v[counter];
          label[u] = iter;
          wl1[i] = u;
          counter++;
        }
      }
      wl2_v.clear();
      if (size1 == 0) break;
      #pragma omp barrier
    } while (true);
    g_level = std::max(g_level, level);
    sum += l_sum;
  }
  dist[v] = g_level - 1;
  return g_level - 1;
}


static int winnow(const int v, int d, const int diameter)
{
  const int max_level = diameter / 2;

  // level-by-level BFS
  iter++;
  label[v] = iter;
  wl1[0] = v;
  int size1 = 1;
  int level = 0;
  int count = 1;
  wl_index[0] = 0;
  #pragma omp parallel default(none) private(level) shared(wl1, g, label, dist, iter, size1, count, wl_index, max_level, diameter) firstprivate(d)
  {
    level = 0;
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    std::vector<int> wl2_v;
    do {
      long start = tid * static_cast<long>(size1) / num_threads;
      long end = (tid + 1) * static_cast<long>(size1) / num_threads;
      level++;
      d++;
      if ((level > max_level) && (d > diameter)) break;
      for (int i = start; i < end; i++) {
        const int u = wl1[i];
        for (int j = g.nindex[u]; j < g.nindex[u + 1]; j++) {
          const int n = g.nlist[j];
          int current;
          #pragma omp atomic capture
          {
            current = label[n];
            label[n] = iter;
          }
          if (current != iter) {
            wl2_v.push_back(n);
            dist[n] = std::min(dist[n], (d <= diameter) ? d : (INT_MAX - 1));
          }
        }
      }
      wl_index[tid + 1] = wl2_v.size();
      #pragma omp barrier
      #pragma omp single
      {
        for (int i = 2; i <= num_threads; i++) {
          count += wl_index[i];
          wl_index[i] += wl_index[i - 1];
        }
        size1 = wl_index[num_threads];
      }
      int counter = 0;
      for (int i = wl_index[tid]; i < wl_index[tid + 1]; i++) {
        wl1[i] = wl2_v[counter];
        counter++;
      }
      wl2_v.clear();
      if (size1 == 0) break;
      #pragma omp barrier
    } while (true);
  }
  return count;
}


static void eliminate(int size1, int d, const int diameter)
{
  // level-by-level BFS
  while (d < diameter) {
    d++;
    int size2 = 0;
    for (int i = 0; i < size1; i++) {
      const int u = wl1[i];
      for (int j = g.nindex[u]; j < g.nindex[u + 1]; j++) {
        const int n = g.nlist[j];
        if (dist[n] > d) {
          dist[n] = d;
          if (d < diameter) {
            wl2[size2] = n;
            size2++;
          }
        }
      }
    }
    if (size2 == 0) break;
    std::swap(wl1, wl2);
    size1 = size2;
  }
}


static void remove_chains(int& diameter)
{
  for (int u = 0; u < g.nodes; u++) {
    if (dist[u] == INT_MAX) {
      const int deg = g.nindex[u + 1] - g.nindex[u];
      if (deg == 1) {  // start of a chain
        int old = dist[u];
        int v = g.nlist[g.nindex[u]];
        int from = u;
        int len = 0;
        do {
          len++;
          const int deg_v = g.nindex[v + 1] - g.nindex[v];
          if (deg_v == 1) {  // end of a linked list
            old = len;
            dist[v] = len;
            diameter = std::max(diameter, len);
            break;
          }
          dist[v] = std::min(dist[v], INT_MAX - 1 - len);  // eliminate
          if (deg_v > 2) {  // end of chain
            wl1[0] = v;
            eliminate(1, INT_MAX - 1 - len, INT_MAX - 1);
            break;
          }
          const int n = g.nlist[g.nindex[v]];
          if (n == from) {  // find next vertex in chain
            from = v;
            v = g.nlist[g.nindex[v] + 1];
          } else {
            from = v;
            v = n;
          }
        } while (true);
        dist[u] = old;
      }
    }
  }
}


int main(int argc, char* argv [])
{
  printf("F-Diam v01 (CPU-parallel graph-diameter computation)\n");
  printf("Copyright 2025 Texas State University\n\n");

  if (argc != 2) {printf("USAGE: %s input_file\n", argv[0]); exit(-1);}

  // read input
  g = readECLgraph(argv[1]);
  printf("input: %s\n", argv[1]);
  printf("nodes: %d\n", g.nodes);
  printf("edges: %d (%d)\n\n", g.edges / 2, g.edges);

  // allocate arrays
  label = new int [g.nodes];
  dist = new int [g.nodes + 1];
  wl1 = new int [g.nodes];
  wl2 = new int [g.nodes];
  #pragma omp parallel
  {
    #pragma omp single
    {
      const int num_threads = omp_get_num_threads();
      wl_index = new int [num_threads + 1];
    }
  }

  // initialize and find highest-degree vertex
  int max_u = 0, max_d = 0;
  for (int u = 0; u < g.nodes; u++) {
    const int deg = g.nindex[u + 1] - g.nindex[u];
    label[u] = -1;
    dist[u] = (deg == 0) ? 0 : INT_MAX;
    if (max_d < deg) {
      max_d = deg;
      max_u = u;
    }
  }
  dist[g.nodes] = INT_MAX;  // sentinel

  CPUTimer timer;
  timer.start();

  // determine initial diameter
  int sum, dummy;
  const int diam = distance(max_u, sum, 0);
  if (sum != g.nodes) printf("graph is disconnected: diameter is infinite\ncomputing maximum diameter of connected components\n");
  const int next = wl1[0];  // a node with maximal distance from max_u
  int diameter = distance(next, dummy, diam);
  printf("initial diameter: %d\n", diameter);
  fflush(stdout);

  // clique detection
  if ((diameter == 1) && (max_d + 1 == sum)) {
    int j = g.nindex[max_u];
    while (j < g.nindex[max_u + 1]) {
      const int n = g.nlist[j];
      const int deg = g.nindex[n + 1] - g.nindex[n];
      if (deg != max_d) break;
      j++;
    }
    if (j == g.nindex[max_u + 1]) {
      for (int j = g.nindex[max_u]; j < g.nindex[max_u + 1]; j++) {
        const int n = g.nlist[j];
        dist[n] = 1;
      }
    }
  }

  // remove "uninteresting" vertices
  const int out = winnow(max_u, diam, diameter);
  printf("vertices reachable from u: %d\n", out);
  fflush(stdout);

  // remove "chain" vertices
  remove_chains(diameter);

  // process vertices
  int v = 0;
  int loop = 0;
  while (true) {
    loop++;
    // determine next starting vertex
    while (dist[v] != INT_MAX) v++;
    if (v == g.nodes) break;

    // compute eccentricity of v
    const int d = distance(v, dummy, diameter);
    if (d > diameter) {
      // new highest bound found
      printf("new bound: %d\n", d);
      fflush(stdout);

      // remove additional "uninteresting" vertices
      if ((d - 1 > diameter) || ((d % 2) == 0)) {
        const int out = winnow(max_u, diam, d);
        printf("vertices reachable from u: %d\n", out);
        fflush(stdout);
      }

      // eliminate additional neighbors
      int s = 0;
      for (int u = 0; u < g.nodes; u++) {
        if (dist[u] == diameter) wl1[s++] = u;
      }
      eliminate(s, diameter, d);
      diameter = d;
    } else {
      // eliminate neighbors
      wl1[0] = v;
      eliminate(1, d, diameter);
    }
  }
  printf("iterations: %d\n", iter);
  printf("final diameter: %d\n\n", diameter);
  printf("runtime: %.4f s\n", timer.elapsed());

  // clean up
  freeECLgraph(g);
  delete [] label;
  delete [] dist;
  delete [] wl1;
  delete [] wl2;
  delete [] wl_index;
  return 0;
}
