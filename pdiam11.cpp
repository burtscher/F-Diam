#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <climits>
#include <set>
#include <random>
#include "ECLgraph.h"
#include <omp.h>
#include <sys/time.h>

static ECLgraph g;
static int iter = 0;
static int* label;
static int* dist;
static int* wl1;
static int* wl2;
static int* wl_index;


/*********************************************************************************/

struct CPUTimer
{
  timeval beg, end;
  CPUTimer() {}
  ~CPUTimer() {}
  void start() {gettimeofday(&beg, NULL);}
  double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};

static void output_graph(const char* const label, const int num, const int mark = INT_MAX)
{
  // output dot file
  char name [256];
  sprintf(name, "output/diam%d.dot", num);
  FILE* f = fopen(name, "w");
  fprintf(f, "graph {\n");
  fprintf(f, "  label = \"%s\";\n", label);
  fprintf(f, "  labelloc = \"t\";\n");
  for (int i = 0; i < g.nodes; i++) {
    if (i == mark) {
      fprintf(f, "  %d [style=filled, fillcolor=\"#ff0000\", label=\"%d\"]\n", i, i);  // red
    } else if (dist[i] != INT_MAX) {
      if (dist[i] > INT_MAX / 2) {
        fprintf(f, "  %d [style=filled, fillcolor=\"#9090ff\", label=\"%d\"]\n", i, i);  // blue
      } else {
        fprintf(f, "  %d [style=filled, fillcolor=\"#00ff00\", label=\"%d\"]\n", i, i);  // green
      }
    } else {
      fprintf(f, "  %d [style=filled, fillcolor=\"#ffffff\", label=\"%d\"]\n", i, i);  // white
    }
  }
  for (int i = 0; i < g.nodes; i++) {
    for (int j = g.nindex[i]; j < g.nindex[i + 1]; j++) {
      const int n = g.nlist[j];
      if (i < n) {
        fprintf(f, "  %d -- %d\n", i, n);
      }
    }
  }
  fprintf(f, "}\n");
  fclose(f);
}


static void output_final(const char* const label)
{
  // output dot file
  FILE* f = fopen("output/final.dot", "w");
  fprintf(f, "graph {\n");
  fprintf(f, "  label = \"%s\";\n", label);
  fprintf(f, "  labelloc = \"t\";\n");
  for (int i = 0; i < g.nodes; i++) {
    fprintf(f, "  %d [style=filled, fillcolor=\"#d0d0d0\", label=\"%d\"]\n", i, dist[i]);
  }
  for (int i = 0; i < g.nodes; i++) {
    for (int j = g.nindex[i]; j < g.nindex[i + 1]; j++) {
      const int n = g.nlist[j];
      if (i < n) {
        fprintf(f, "  %d -- %d\n", i, n);
      }
    }
  }
  fprintf(f, "}\n");
  fclose(f);
}


/*********************************************************************************/


static int distance(const int v, int& sum, int diameter)
{
  // level-by-level BFS
  iter++;
  label[v] = iter;
  wl1[0] = v;
  int g_level = 0;
  int size1 = g.nindex[v + 1] - g.nindex[v];
  sum = size1 + 1;
  wl_index[0] = 0;
  #pragma omp parallel default(none) shared(wl1, g, label, iter, size1, wl_index, v) reduction(max:g_level) reduction(+:sum)
  {
    int level = 0;
    int l_sum = 0;
    int curVal;
    int counter = 0;
    int tid = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    std::vector<int> wl2_v;
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
    float ratio;
    level++;
    //for loop through first adjacency list
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
        //wl pre-work
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
          for (int i = 1; i <= num_threads; i++) {
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
            for (int i = 1; i <= num_threads; i++) {
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
        for (int i = 1; i <= num_threads; i++) {
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
  int cnt = 0;  //MB: stats only
  int max = 0;  //MB: stats only
  for (int u = 0; u < g.nodes; u++) {
    if (dist[u] == INT_MAX) {
      const int deg = g.nindex[u + 1] - g.nindex[u];
      if (deg == 1) {  // start of a chain
        cnt++;  //MB: stats only
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
        max = std::max(max, len);  //MB: stats only
      }
    }
  }
  printf("found %d chains (max = %d)\n", cnt, max);  //MB: stats only
  fflush(stdout);  //MB: stats only
}


int main(int argc, char* argv [])
{
  printf("CPU parallel graph diameter v01 (%s)\n", __FILE__);
  printf("Copyright 2024 Texas State University\n\n");

  if (argc < 2) {printf("USAGE: %s input_file\n", argv[0]); exit(-1);}

  if (argc == 2) {

    // read input
    g = readECLgraph(argv[1]);
    printf("input: %s\n", argv[1]);
    printf("nodes: %d\n", g.nodes);
    printf("edges: %d (%d)\n", g.edges / 2, g.edges);
  } else {
    // process command line
    if ((argc < 3) || ((argc % 2) != 1)) {fprintf(stderr, "USAGE: %s number_of_vertices number_of_edges [extra_edge_pair ...]\n", argv[0]); exit(-1);}
    const int n = atoi(argv[1]);
    if (n < 2) {fprintf(stderr, "ERROR: number_of_vertices must be at least 2\n"); exit(-1);}
    int m = atoi(argv[2]);
    if ((m < 0) || (m > n * (n - 1) / 2)) {fprintf(stderr, "ERROR: number_of_edges must be between 0 and %d\n", n * (n - 1) / 2); exit(-1);}

    // generate random graph
    std::set<int>* adj = new std::set<int> [n];
    std::uniform_int_distribution<int> rndval(0, n - 1);
    std::mt19937 mt(23);
    for (int j = 0; j < m; j++) {
      int src, dst;
      do {
        do {
          src = rndval(mt);
          dst = rndval(mt);
        } while (src == dst);
      } while (adj[dst].find(src) != adj[dst].end());
      adj[dst].insert(src);
      adj[src].insert(dst);
    }
    m *= 2;

    // process extra edges
    int p = 3;
    while (p + 1 < argc) {
      const int s = atoi(argv[p]);
      const int d = atoi(argv[p + 1]);
      if ((s < 0) || (s >= n)) {fprintf(stderr, "ERROR: source must be between 0 and %d\n", n - 1); exit(-1);}
      if ((d < 0) || (d >= n)) {fprintf(stderr, "ERROR: destination must be between 0 and %d\n", n - 1); exit(-1);}
      if (adj[s].find(d) != adj[s].end()) {
        adj[s].erase(d);
        adj[d].erase(s);
      } else {
        adj[s].insert(d);
        adj[d].insert(s);
      }
      m += 2;
      p += 2;
    }
    printf("%d nodes and %d edges\n", n, m);

    g.nodes = n;
    g.edges = m;
    g.nindex = new int [n + 1];
    g.nlist = new int [m];
    g.eweight = NULL;

    int pos = 0;
    for (int i = 0; i < n; i++) {
      g.nindex[i] = pos;
      for (auto v: adj[i]) {
        g.nlist[pos] = v;
        pos++;
      }
    }
    g.nindex[n] = pos;

    delete [] adj;
  }

  // allocate arrays
  label = new int [g.nodes];
  dist = new int [g.nodes + 1];
  wl1 = new int [g.nodes];
  wl2 = new int [g.nodes];
  #pragma omp parallel
  {
    #pragma single
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

  if (g.nodes < 40) output_graph("Highest degree node", 0, max_u);
  CPUTimer timer;
  timer.start();
  // determine initial diameter
  int sum, dummy;
  const int diam = distance(max_u, sum, 0);
  if (sum != g.nodes) printf("graph is disconnected: diameter is infinite\n  computing maximum diameter of connected components %d %d\n", sum, g.nodes);
  printf("max_u %d: %d\n", max_u, diam);
  const int next = wl1[0];  // a node with maximal distance from max_u
  int diameter = distance(next, dummy, diam);
  printf("initial diameter: %d\n", diameter);
  fflush(stdout);
  if (g.nodes < 40) output_graph("Farthest node", 1, next);
  if (g.nodes < 40) output_graph("Next farthest node", 2, wl1[0]);

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

  if (g.nodes < 40) output_graph("After clique detection", 3);

  // remove "uninteresting" vertices
  const int out = winnow(max_u, diam, diameter);
  printf("d/2 reachable from max_u: %d\n", out);
  fflush(stdout);

  if (g.nodes < 40) output_graph("After removing uninteresting nodes", 4);

  // remove "chain" vertices
  remove_chains(diameter);

  if (g.nodes < 40) output_graph("After removing chains", 5);

  // process vertices
  int v = 0;
  int loop = 0;
  while (true) {
    loop++;
    // determine next starting vertex
    while (dist[v] != INT_MAX) v++;
    if (v == g.nodes) break;

    if (g.nodes < 40) output_graph("Current starting vertex", loop * 10 + 0, v);

    // compute eccentricity of v
    const int d = distance(v, dummy, diameter);

    if (d > diameter) {
      // new highest bound found
      printf("diam = %d (node = %d to %d, iter = %d)\n", d, v, wl1[0], iter);
      fflush(stdout);

      // remove additional "uninteresting" vertices
      if ((d - 1 > diameter) || ((d % 2) == 0)) {
        const int out = winnow(max_u, diam, d);
        printf("d/2 reachable from max_u: %d\n", out);
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

    if (g.nodes < 40) output_graph("After eliminating nodes", loop * 10 + 1);
  }
  printf("iterations: %d\n", iter);
  printf("final diameter: %d\n\n", diameter);
  printf("runtime: %.4f s\n", timer.elapsed());
  if (g.nodes < 40) {
    for (int v = 0; v < g.nodes; v++) {
      distance(v, dummy, diameter);  // actually compute all eccentricities
    }
    output_final("Eccentricities");
  }
  // clean up
  freeECLgraph(g);
  delete [] label;
  delete [] dist;
  delete [] wl1;
  delete [] wl2;
  delete [] wl_index;
  return 0;
}


/*
rm -f output/*; ./dia44 9 13
rm -f output/*; ./dia44 12 16 3 11
rm -f output/*; ./dia44 14 16 4 7
rm -f output/*; ./dia44 24 33
rm -f output/*; ./dia44 10 20 8 5 8 1 2 0 2 6 2 7 2 9
for f in output/*.dot; do dot -Tpng $f > $f.png; done
*/


