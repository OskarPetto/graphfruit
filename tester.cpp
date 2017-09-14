#include "graph.hpp"
#include "digraph.hpp"

using namespace graphfruit;

int main(int argc, char const *argv[]) {

  graph<> g;

  g.add_edge(0, 1, 4.0);
  g.add_edge(1, 2, 8.0);
  g.add_edge(2, 3, 7.0);
  g.add_edge(3, 4, 9.0);
  g.add_edge(4, 5, 10.0);
  g.add_edge(5, 6, 2.0);
  g.add_edge(6, 7, 1.0);
  g.add_edge(7, 0, 8.0);
  g.add_edge(7, 1, 11.0);
  g.add_edge(7, 8, 7.0);
  g.add_edge(8, 2, 2.0);
  g.add_edge(8, 6, 6.0);
  g.add_edge(2, 5, 4.0);
  g.add_edge(3, 5, 14.0);

  digraph<> d1(6);

  d1.add_edge(0, 1, 5);
  d1.add_edge(0, 2, 3);
  d1.add_edge(1, 2, 2);
  d1.add_edge(1, 3, 6);
  d1.add_edge(2, 3, 7);
  d1.add_edge(2, 4, 4);
  d1.add_edge(2, 5, 2);
  d1.add_edge(3, 4, -1);
  d1.add_edge(3, 5, 1);
  d1.add_edge(4, 5, -2);

  digraph<> d2(5);

  d2.add_edge(0, 1, -1);
  d2.add_edge(0, 2, 4);
  d2.add_edge(1, 2, 3);
  d2.add_edge(1, 3, 2);
  d2.add_edge(1, 4, 2);
  d2.add_edge(3, 2, 5);
  d2.add_edge(3, 1, 1);
  d2.add_edge(4, 3, -3);

  std::cout << " -------Graph------- " << '\n';

  std::cout << g.degree(4) << '\n';
  std::cout << g.edge_weight(0, 1) << '\n';

  std::cout << " -------Digraph------- " << '\n';

  std::cout << d1.in_degree(4) << '\n';
  std::cout << d1.out_degree(4) << '\n';
  std::cout << d1.edge_weight(2, 3) << '\n';

  std::cout << " -------Breadth-first search------- " << '\n';

  std::vector<bool> visited = breadth_first_search(d1, 2);

  for (bool i : visited) {
    std::cout << i << '\n';
  }

  std::cout << " -------Shortest path------- " << '\n';

  size_t start = 0;
  size_t end = 4;

  std::vector<size_t> path = dijkstra_shortest_path(g, start, end);

  for (size_t i : path) {
    std::cout << i << '\n';
  }

  std::cout << " -------Shortest path (negative)------- " << '\n';

  start = 0;
  end = 3;

  path = bellman_ford_shortest_path(d2, start, end);

  for (size_t i : path) {
    std::cout << i << '\n';
  }

  std::cout << " -------All shortest paths------- " << '\n';

  std::vector<std::vector<ssize_t> > all_paths = johnson_all_shortest_paths(g);

  for (ssize_t i = 0; i < all_paths.size(); i++) {
    std::cout << i << ": ";
    for (ssize_t j : all_paths[i]) {
      std::cout << j << " ";
    }
    std::cout << '\n';
  }

  std::cout << " -------Is cyclic------- " << '\n';

  bool cycle = is_cyclic(d1);
  std::cout << cycle << '\n';

  std::cout << " -------Minimum spanning tree------- " << '\n';

  std::vector<std::pair<size_t, size_t> > kmst = kruskal_minimum_spanning_tree(g);

  for (std::pair<size_t, size_t> i : kmst) {
    std::cout << "(" << i.first << ", " << i.second << ")" << '\n';
  }

  std::cout << " -------Topological sorting------- " << '\n';

  std::vector<size_t> sorting = khan_topological_sort(d1);

  for (size_t i : sorting) {
    std::cout << i << '\n';
  }

  std::cout << " -------Longest path in DAG------- " << '\n';

  path = DAG_longest_path(d1, 1, 5);

  for (size_t i : path) {
    std::cout << i << '\n';
  }

  return 0;
}
