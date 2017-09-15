/*
 * A class for undirected graphs.
 * @version 15.09.2017
    added k-core and changed formating of operator <<
 * @version 14.09.2017
 *  Johnson's algorithm added
 * @version 13.09.2017
 *  first version
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include "base_graph.hpp"
#include "union_find.hpp"

namespace graphfruit {

  /*
   * An object of this class represents an undirected graph.
   * V is the type of data that can be stored in a vertex.
   */
  template <class V = void*>
  class graph : public base_graph<V> {

    typedef typename graph<V>::edge edge;
    typedef typename graph<V>::vertex vertex;

  public:

    /*
     * Constructors
     */
    explicit graph(size_t n = 0);
    graph(const graph<V>& other) = default;
    graph(graph<V>&& other) noexcept = default;
    graph& operator=(const graph<V>& other) = default;
    graph& operator=(graph<V>&& other) noexcept = default;
    ~graph() = default;

    /*
     * Makes a digraph printable to an output stream.
     */
    template <class V1>
    friend std::ostream& operator<<(std::ostream& out, const graph<V1>& g);

    /*
     * Adds an undirected edge between the source vertex and the target vertex
     * to the graph and adds the vertices if they don't exist in the graph.
     * Complexity: O(1) amortized
     */
    void add_edge(size_t source_vertex, size_t target_vertex, double edge_weight = 1.0);

    /*
     * Returns the degree of the vertex. Returns 0 if the vertex does not exist.
     * Complexity: O(1)
     */
    size_t degree(size_t u) const;

    /*
     * Returns a vector of vectex pairs representing all edges.
     * Complexity: O(E)
     */
    std::vector<std::pair<size_t, size_t> > edges() const;

    /*
     * Returns the neighbour vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(E)
     */
    std::vector<size_t> neighbours(size_t u) const;

    /*
     * Returns the number of undirected edges in the graph.
     * Complexity: O(1)
     */
    size_t number_of_edges() const {return this->edge_list.size() / 2;}

    /*
     * Removes all edges between the source vertex and the target vertex from
     * the graph. Does nothing if there are no undirected edges.
     * Complexity: O(E^2)
     */
    void remove_edges(size_t source_vertex, size_t target_vertex);

    /*
     * Uses Johnsons's algorithm to calculate the shortest paths between all
     * pairs of vertices. Returns a 2D vector of predecessors in these shortest
     * paths. If there exists no path -1 is the predecessor. Returns an empty
     * 2D vector if the graph contains a negative cycle.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::vector<ssize_t> > johnson_all_shortest_paths(const graph<V1>& g);

    /*
     * Returns the k-core of the graph. A k-core is a maximal subgraph in which
     * all vertices have a degree of k or more.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend graph<V1> k_core(const graph<V1>& g, size_t k);

    /*
     * Uses Kruskal's algorithm to find the edges of a minimum spanning tree of
     * the graph. Returns the edges a vector of vertex pairs. If the graph is
     * not connnected the result is not defined.
     * Complexity: O(E * log(V))
     */
    template <class V1>
    friend std::vector<std::pair<size_t, size_t> > kruskal_minimum_spanning_tree(const graph<V1>& g);

    /*
     * Uses Prims's algorithm to find the edges of a minimum spanning tree of
     * the graph. Returns the edges a vector of vertex pairs. If the graph is
     * not connnected the result is not defined.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::pair<size_t, size_t> > prim_minimum_spanning_tree(const graph<V1>& g);

  protected:

    // for sort in kruskal_minimum_spanning_tree
    struct edge_comp {
      bool operator() (edge* a, edge*b) {
        return a->edge_weight <= b->edge_weight;
      }
    };

    bool DFS_degree(std::vector<bool>& visited, std::size_t u, std::vector<std::size_t>& degree, std::size_t k) const;

  };

  /*
   * graph member implementations
   */
  template <class V>
  graph<V>::graph(size_t n)
    : base_graph<V>(n) {}

  template <class V>
  std::ostream& operator<<(std::ostream& out, const graph<V>& g) {
    out << "Graph: V=" << g.number_of_vertices();
    out << ", E=" << g.number_of_edges();
    // for (size_t i = 0; i < g.number_of_edges(); i++) {
    //   out << std::endl;
    //   out << " ";
    //   out << "(" << g.edge_list[2 * i]->source_vertex->vertex_index;
    //   out << ", " << g.edge_list[2 * i]->target_vertex->vertex_index;
    //   out << ", " << g.edge_list[2 * i]->edge_weight;
    //   out << ")";
    // }
    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      out << std::endl;
      out << " ";
      out << "[" << g.vertex_list[i]->vertex_index;
      out << "]";
      for (typename graph<V>::edge* e : g.vertex_list[i]->outgoing_edge_list) {
        out << " - " << e->target_vertex->vertex_index;
      }
    }
    return out;
  }

  template <class V>
  void graph<V>::add_edge(size_t source_vertex, size_t target_vertex, double edge_weight) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      size_t max = source_vertex > target_vertex ? source_vertex + 1: target_vertex + 1;
      size_t i = this->number_of_vertices();
      this->vertex_list.resize(max);
      for (;i < max; i++) {
        this->vertex_list[i] = new vertex(i);
      }
    }
    edge* e1 = new edge(this->vertex_list[source_vertex], this->vertex_list[target_vertex], edge_weight);
    edge* e2 = new edge(this->vertex_list[target_vertex], this->vertex_list[source_vertex], edge_weight);
    this->edge_list.push_back(e1);
    this->edge_list.push_back(e2);
    this->vertex_list[source_vertex]->outgoing_edge_list.push_back(e1);
    this->vertex_list[target_vertex]->outgoing_edge_list.push_back(e2);
  }

  template <class V>
  size_t graph<V>::degree(size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    return this->vertex_list[u]->outgoing_edge_list.size();
  }


  template <class V>
  std::vector<std::pair<size_t, size_t> > graph<V>::edges() const {
    std::vector<std::pair<size_t, size_t> > v(number_of_edges());
    for (size_t i = 0; i < number_of_edges(); i++) {
      std::pair<size_t, size_t> a;
      a.first = this->edge_list[i * 2]->source_vertex->vertex_index;
      a.second = this->edge_list[i * 2]->target_vertex->vertex_index;
      v[i] = a;
    }
    return v;
  }

  template <class V>
  std::vector<size_t> graph<V>::neighbours(size_t u) const {
    std::vector<size_t> v;
    if (!this->contains_vertex(u)) {
      return v;
    }
    for (edge* e : this->vertex_list[u]->outgoing_edge_list) {
      v.push_back(e->target_vertex->vertex_index);
    }
    return v;
  }

  template <class V>
  void graph<V>::remove_edges(size_t source_vertex, size_t target_vertex) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      return;
    }
    typename std::vector<edge*>::iterator i = this->edge_list.begin();
    while (i != this->edge_list.end()) {
      if (((*i)->source_vertex == this->vertex_list[source_vertex] &&
           (*i)->target_vertex == this->vertex_list[target_vertex]    )
          || ((*i)->source_vertex == this->vertex_list[target_vertex] &&
              (*i)->target_vertex == this->vertex_list[source_vertex]    )) {
        delete *i;
        i = this->edge_list.erase(i);
      } else {
        ++i;
      }
    }
    typename std::vector<edge*>::iterator j = this->vertex_list[source_vertex]->outgoing_edge_list.begin();
    while (j != this->vertex_list[source_vertex]->outgoing_edge_list.end()) {
      if ((*j)->target_vertex == this->vertex_list[target_vertex]) {
        j = this->vertex_list[source_vertex]->outgoing_edge_list.erase(j);
      } else {
        ++j;
      }
    }
    typename std::vector<edge*>::iterator k = this->vertex_list[target_vertex]->outgoing_edge_list.begin();
    while (k != this->vertex_list[target_vertex]->outgoing_edge_list.end()) {
      if ((*k)->target_vertex == this->vertex_list[source_vertex]) {
        k = this->vertex_list[target_vertex]->outgoing_edge_list.erase(k);
      } else {
        ++k;
      }
    }
  }

  template <class V>
  std::vector<std::vector<ssize_t> > johnson_all_shortest_paths(const graph<V>& g) {
    graph<V> g1(g);
    size_t u = g1.number_of_vertices();
    g1.add_vertex();
    for (size_t v = 0; v < u; v++) {
      g1.add_edge(u, v, 0.0);
    }
    std::vector<double> bf_distance = g1.bellman_ford_distance(u);
    if (bf_distance.empty()) {
      std::vector<std::vector<ssize_t> > empty;
      return empty;
    }
    g1 = g;
    for (typename base_graph<V>::edge* e : g1.edge_list) {
      size_t u = e->source_vertex->vertex_index;
      size_t v = e->target_vertex->vertex_index;
      e->edge_weight += bf_distance[u] - bf_distance[v];
    }
    std::vector<std::vector<ssize_t> > previous(g1.number_of_vertices());
    for (size_t i = 0; i < g1.number_of_vertices(); i++) {
      previous[i] = dijkstra_shortest_path(g1, i);
    }
    return previous;
  }

  template <class V>
  graph<V> k_core(const graph<V>& g, size_t k) {
    std::vector<bool> visited(g.number_of_vertices());
    std::vector<size_t> degree(g.number_of_vertices());

    size_t min_degree = -1;
    size_t start_vertex;

    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      degree[i] = g.vertex_list[i]->outgoing_edge_list.size();
      if (degree[i] < min_degree) {
        min_degree = degree[i];
        start_vertex = i;
      }
    }
    g.DFS_degree(visited, start_vertex, degree, k);
    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      if (!visited[i]) {
        g.DFS_degree(visited, i, degree, k);
      }
    }
    std::vector<ssize_t> pos(g.number_of_vertices(), -1);
    size_t i = 0;
    graph<V> g1;
    for (size_t u = 0; u < g.number_of_vertices(); u++) {
      if (degree[u] >= k) {
        g1.add_vertex(g.vertex_list[u]->vertex_data);
        pos[u] = i;
        i++;
      }
    }
    for (typename graph<V>::edge* e : g.edge_list) {
      size_t u = e->source_vertex->vertex_index;
      size_t v = e->target_vertex->vertex_index;
      if (degree[u] >= k && degree[v] >= k && u > v) {
        g1.add_edge(pos[u], pos[v], e->edge_weight);
      }
    }
    return g1;
  }

  template <class V>
  bool graph<V>::DFS_degree(std::vector<bool>& visited, size_t u, std::vector<size_t>& degree, size_t k) const {
    visited[u] = true;

    for (edge* e : this->vertex_list[u]->outgoing_edge_list) {
      size_t v = e->target_vertex->vertex_index;
      if (degree[u] < k) {
        degree[v]--;
      }
      if (!visited[v]) {
        if (DFS_degree(visited, v, degree, k)) {
          degree[u]--;
        }
      }
    }
    return degree[u] < k;
  }

  template <class V>
  std::vector<std::pair<size_t, size_t> > kruskal_minimum_spanning_tree(const graph<V>& g) {

    std::vector<std::pair<size_t, size_t> > result(g.number_of_vertices() - 1);
    std::vector<typename graph<V>::edge*> sorted_edges(g.number_of_edges());
    for (size_t i = 0; i < g.number_of_edges(); i++) {
      sorted_edges[i] = g.edge_list[2 * i];
    }
    std::sort(sorted_edges.begin(), sorted_edges.end(), typename graph<V>::edge_comp());
    union_find uf(g.number_of_vertices());
    size_t i = 0;
    for (typename graph<V>::edge* e : sorted_edges) {
      if (i == g.number_of_vertices() - 1) {
        return result;
      }
      size_t u_root = uf.set_find(e->source_vertex->vertex_index);
      size_t v_root = uf.set_find(e->target_vertex->vertex_index);
      if (u_root != v_root) {
        result[i] = std::make_pair(e->source_vertex->vertex_index, e->target_vertex->vertex_index);
        i++;
        uf.set_union(u_root, v_root);
      }
    }
    return result;
  }

  template <class V>
  std::vector<std::pair<size_t, size_t> > prim_minimum_spanning_tree(const graph<V>& g) {

    std::vector<std::pair<size_t, size_t> > result(g.number_of_vertices() - 1);
    std::vector<double> distance(g.number_of_vertices(), std::numeric_limits<double>::max());
    std::vector<fibonacci_node<std::pair<size_t, double> >*> fib_nodes(g.number_of_vertices());
    fibonacci_heap<std::pair<size_t, double>, typename graph<V>::fib_comp> min_heap;

    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      std::pair<size_t, double> a(i, std::numeric_limits<double>::max());
      fib_nodes[i] = min_heap.push(a);
    }
    distance[0] = 0.0;
    min_heap.decrease_key(fib_nodes[0], std::make_pair(0, 0.0));
    while (!min_heap.empty()) {
      std::pair<size_t, double> u = min_heap.top();
      fib_nodes[u.first] = nullptr;
      min_heap.pop();
      for (typename graph<V>::edge* e : g.vertex_list[u.first]->outgoing_edge_list) {
        size_t v = e->target_vertex->vertex_index;
        if (fib_nodes[v] && distance[v] > e->edge_weight) {
          distance[v] = e->edge_weight;
          min_heap.decrease_key(fib_nodes[v], std::make_pair(v, distance[v]));
          result[v-1] = std::make_pair(u.first, v);
        }
      }
    }
    return result;
  }

}

#endif //GRAPH_HPP
