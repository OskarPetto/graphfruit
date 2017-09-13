#include "graph.hpp"
#include "digraph.hpp"
#include <iostream>
#include <type_traits>

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

  digraph<> d(6);

  d.add_edge(5, 2);
  d.add_edge(5, 0);
  d.add_edge(4, 0);
  d.add_edge(4, 1);
  d.add_edge(2, 3);
  d.add_edge(3, 1);

  std::cout << " -------Breadth First Search------- " << '\n';

  std::vector<bool> visited = breadth_first_search(d, 2);

  for (bool i : visited) {
    std::cout << i << '\n';
  }

  std::cout << " -------Shortest Path------- " << '\n';

  std::size_t start = 0;
  std::size_t end = 4;

  std::vector<std::size_t> path = dijkstra_shortest_path(g, start, end);

  for (std::size_t i : path) {
    std::cout << i << '\n';
  }

  std::cout << " -------Is Cyclic------- " << '\n';

  bool cycle = is_cyclic(d);
  std::cout << cycle << '\n';

  std::cout << " -------Minimum Spanning Tree------- " << '\n';

  std::vector<std::pair<std::size_t, std::size_t>> kmst = kruskal_minimum_spanning_tree(g);
  std::vector<std::pair<std::size_t, std::size_t>> pmst = prim_minimum_spanning_tree(g);

  for (std::size_t i = 0; i < kmst.size(); i++) {
    std::cout << "(" << kmst[i].first << ", " << kmst[i].second << ")";
    std::cout << " (" << pmst[i].first << ", " << pmst[i].second << ")" << '\n';
  }

  std::cout << " -------Topological Sort------- " << '\n';

  std::vector<std::size_t> sorting = khan_topological_sort(d);

  for (std::size_t i : sorting) {
    std::cout << i << '\n';
  }

  return 0;
}
