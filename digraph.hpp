/*
 * A class for directed graphs.
 * @version 21.09.2017
 *  changed vector.resize to vector.reserve
 * @version 18.09.2017
 *  subgraph and subgraphs added
 * @version 16.09.2017
 *  transpose added
 * @version 15.09.2017
 *  changed formating of operator <<
 * @version 14.09.2017
 *  Johnson's algorithm added
 * @version 13.09.2017
 *  first version
 */

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include "base_graph.hpp"

namespace graphfruit {

  /*
   * An object of this class represents a directed graph.
   * V is the type of data that can be stored in a vertex.
   */
  template <class V = void*>
  class digraph : public base_graph<V> {

    typedef typename digraph<V>::edge edge;
    typedef typename digraph<V>::vertex vertex;

  public:

    /*
     * Constructors
     */
    explicit digraph(size_t n = 0);
    digraph(const digraph<V>& other) = default;
    digraph(digraph<V>&& other) noexcept = default;
    digraph& operator=(const digraph<V>& other) = default;
    digraph& operator=(digraph<V>&& other) noexcept = default;
    ~digraph() = default;

    /*
     * Makes a digraph printable to an output stream.
     */
    template <class V1>
    friend std::ostream& operator<<(std::ostream& out, const digraph<V1>& g);

    /*
     * Adds a directed edge between the source vertex and the target vertex to
     * the graph and adds the vertices if they don't exist in the graph.
     * Complexity: O(1) amortized
     */
    void add_edge(size_t source_vertex, size_t target_vertex, double edge_weight = 1.0);

    /*
     * Returns a vector of vectex pairs representing all edges.
     * Complexity: O(E)
     */
    std::vector<std::pair<size_t, size_t> > edges() const;

    /*
     * Returns the indegree of the vertex. Returns 0 if the vertex doesn't
     * exist.
     * Complexity: O(E)
     */
    size_t in_degree(size_t u) const;

    /*
     * Returns the number of directed edges in the graph.
     * Complexity: O(1)
     */
    size_t number_of_edges() const {return this->edge_list.size();}

    /*
     * Returns the outdegree of the vertex. Returns 0 if the vertex doesn't
     * exist.
     * Complexity: O(1)
     */
    size_t out_degree(size_t u) const;

    /*
     * Returns the predecessor vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(E)
     */
    std::vector<size_t> predecessors(size_t u) const;

    /*
     * Removes all edges between the source vertex and the target vertex from
     * the graph. Does nothing if there are no directed edges.
     * Complexity: O(E^2)
     */
    void remove_edges(size_t source_vertex, size_t target_vertex);

    /*
     * Returns the successor vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(V)
     */
    std::vector<size_t> successors(size_t u) const;

    /*
     * Calculates the longest paths between the start vertex to all other
     * vertices in a directed acyclic graph (DAG). Returns a vector of
     * predecessors in these longest paths. If there exists no path -1 is the
     * predecessor. Returns an empty vector if the start vertex is not in the
     * graph or the graph is no DAG.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<ssize_t> DAG_longest_path(const digraph<V1>& g, size_t start_vertex);

    /*
     * Calculates the longest path between the start vertex to all other
     * vertices in a directed acyclic graph (DAG). Returns a vector of vertices
     * in the path in reversed order. Returns an empty vector if the start
     * vertex is not in the graph or the graph is no DAG.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<size_t> DAG_longest_path(const digraph<V1>& g, size_t start_vertex, size_t end_vertex);

    /*
     * Returns true if and only if the directed graph is acyclic.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend bool is_DAG(const digraph<V1>& g);

    /*
     * Uses Johnsons's algorithm to calculate the shortest paths between all
     * pairs of vertices. Returns a 2D vector of predecessors in these shortest
     * paths. If there exists no path -1 is the predecessor. Returns an empty
     * 2D vector if the graph contains a negative cycle.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::vector<ssize_t> > johnson_all_shortest_paths(const digraph<V1>& g);

    /*
     * Uses Khan's algorithm to get a topological sorting of the graph. Returns
     * a vector of ordered vertices. If no such sorting is possible an empty
     * vector is returned.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<size_t> khan_topological_sort(const digraph<V1>& g);

    /*
     * Returns a subgraph from an input directed graph and a bool vector
     * indicating which vertices are in the subgraph. If the length of the
     * vector doesn't equal the number of vertices in the graph, an
     * invalid_argument exception is thrown.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend digraph<V1> subgraph(const digraph<V1>& g, std::vector<bool> contains);

    /*
     * Returns a vector of subgraphs from an input directed graph and a
     * size_t vector indicating which vertices belong in which subgraph. If the
     * length of the vector doesn't equal the number of vertices in the graph,
     * an invalid_argument exception is thrown.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<digraph<V1> > subgraphs(const digraph<V1>& g, std::vector<size_t> components);

    /*
     * Returns a directed graph with all of the edges reversed compared to the
     * original directed graph.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend digraph<V1> transpose(const digraph<V1>& g);

  };

  /*
   * digraph member implementations
   */
  template <class V>
  digraph<V>::digraph(size_t n)
    : base_graph<V>(n) {}

  template <class V>
  std::ostream& operator<<(std::ostream& out, const digraph<V>& g) {
    out << "Digraph: V=" << g.number_of_vertices();
    out << ", E=" << g.number_of_edges();
    // for (typename digraph<V>::edge* e : g.edge_list) {
    //   out << std::endl;
    //   out << " ";
    //   out << "(" << e.source_vertex->vertex_index;
    //   out << ", " << e.target_vertex->vertex_index;
    //   out << ") " << e.edge_weight;
    // }
    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      out << std::endl;
      out << " ";
      out << "[" << g.vertex_list[i]->vertex_index;
      out << "]";
      for (typename digraph<V>::edge* e : g.vertex_list[i]->outgoing_edge_list) {
        out << " - " << e->target_vertex->vertex_index;
      }
    }
    return out;
  }

  template <class V>
  void digraph<V>::add_edge(size_t source_vertex, size_t target_vertex, double edge_weight) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      size_t max = source_vertex > target_vertex ? source_vertex + 1 : target_vertex + 1;
      size_t i = this->number_of_vertices();
      this->vertex_list.reserve(max);
      for (;i < max; i++) {
        this->vertex_list.push_back(new vertex(i));
      }
    }
    edge* e = new edge(this->vertex_list[source_vertex], this->vertex_list[target_vertex], edge_weight);
    this->edge_list.push_back(e);
    this->vertex_list[source_vertex]->outgoing_edge_list.push_back(e); // seg fault
  }


  template <class V>
  std::vector<std::pair<size_t, size_t> > digraph<V>::edges() const {
    std::vector<std::pair<size_t, size_t> > v(number_of_edges());
    for (size_t i = 0; i < number_of_edges(); i++) {
      std::pair<size_t, size_t> a;
      a.first = this->edge_list[i]->source_vertex->vertex_index;
      a.second = this->edge_list[i]->target_vertex->vertex_index;
      v[i] = a;
    }
    return v;
  }

  template <class V>
  size_t digraph<V>::in_degree(size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    size_t d = 0;
    for (edge* e : this->edge_list) {
      if (e->target_vertex->vertex_index == u) {
        d++;
      }
    }
    return d;
  }

  template <class V>
  size_t digraph<V>::out_degree(size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    return this->vertex_list[u]->outgoing_edge_list.size();
  }

  template <class V>
  std::vector<size_t> digraph<V>::predecessors(size_t u) const {
    std::vector<size_t> v;
    if (!this->contains_vertex(u)) {
      return v;
    }
    for (edge* e : this->edge_list) {
      if (e->target_vertex->vertex_index == u) {
        v.push_back(e->source_vertex->vertex_index);
      }
    }
    return v;
  }

  template <class V>
  void digraph<V>::remove_edges(size_t source_vertex, size_t target_vertex) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      return;
    }
    typename std::vector<edge*>::iterator i = this->edge_list.begin();
    while (i != this->edge_list.end()) {
      if ((*i)->source_vertex == this->vertex_list[source_vertex] &&
          (*i)->target_vertex == this->vertex_list[target_vertex]) {
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
  }

  template <class V>
  std::vector<size_t> digraph<V>::successors(size_t u) const {
    std::vector<size_t> v;
    if (!this->contains_vertex(u)) {
      return v;
    }
    for (edge* e : this->vertex_list[u].outgoing_edge_list) {
      v.push_back(e->target_vertex->vertex_index);
    }
    return v;
  }

  template <class V>
  std::vector<ssize_t> DAG_longest_path(const digraph<V>& g, size_t start_vertex) {
    if (!g.contains_vertex(start_vertex)) {
      std::vector<ssize_t> empty;
      return empty;
    }
    std::vector<size_t> top_sort = khan_topological_sort(g);
    if (top_sort.empty()) {
      std::vector<ssize_t> empty;
      return empty;
    }
    std::vector<ssize_t> previous(g.number_of_vertices(), -1);
    std::vector<double> distance(g.number_of_vertices(), std::numeric_limits<double>::min());
    distance[start_vertex] = 0.0;

    for (size_t u : top_sort) {
      if (distance[u] != std::numeric_limits<double>::min()) {
        for (typename digraph<V>::edge* e : g.vertex_list[u]->outgoing_edge_list) {
          size_t v = e->target_vertex->vertex_index;
          if (distance[u] + e->edge_weight > distance[v]) {
            distance[v] = distance[u] + e->edge_weight;
            previous[v] = u;
          }
        }
      }
    }
    return previous;
  }

  template <class V>
  std::vector<size_t> DAG_longest_path(const digraph<V>& g, size_t start_vertex, size_t end_vertex) {
    std::vector<size_t> path;
    if (!g.contains_vertex(start_vertex) || !g.contains_vertex(end_vertex)) {
      return path;
    }
    if (start_vertex == end_vertex) {
      return path;
    }
    std::vector<ssize_t> previous = DAG_longest_path(g, start_vertex);
    if (previous.empty() || previous[end_vertex] == -1) {
      return path;
    }
    size_t i = end_vertex;
    while (i != start_vertex) {
      path.push_back(i);
      i = previous[i];
    }
    path.push_back(start_vertex);
    return path;
  }

  template <class V>
  bool is_DAG(const digraph<V>& g) {
    if (g.empty()) {
      return true;
    }
    std::vector<size_t> top_sort = khan_topological_sort(g);
    return !top_sort.empty();
  }

  template <class V>
  std::vector<std::vector<ssize_t> > johnson_all_shortest_paths(const digraph<V>& g) {
    digraph<V> g1(g);
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
  std::vector<size_t> khan_topological_sort(const digraph<V>& g) {
    std::vector<size_t> sorted_vertices(g.number_of_vertices());
    std::vector<size_t> count(g.number_of_vertices());
    std::queue<size_t> q;
    for (typename digraph<V>::edge* e : g.edge_list) {
      count[e->target_vertex->vertex_index]++;
    }
    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      if (count[i] == 0) {
        q.push(i);
      }
    }
    size_t i = 0;
    while (!q.empty()) {
      size_t u = q.front();
      q.pop();
      sorted_vertices[i] = u;
      i++;
      for (typename digraph<V>::edge* e : g.vertex_list[u]->outgoing_edge_list) {
        size_t v = e->target_vertex->vertex_index;
        count[v]--;
        if (count[v] == 0) {
          q.push(v);
        }
      }
    }
    if (i != g.number_of_vertices()) {
      std::vector<size_t> empty;
      return empty;
    }
    return sorted_vertices;;
  }

  template <class V>
  digraph<V> subgraph(const digraph<V>& g, std::vector<bool> contains) {
    if (g.number_of_vertices() != contains.size()) {
      throw std::invalid_argument("subgraph::vector size differs from number of vertices");
    }
    digraph<V> g1;
    std::vector<size_t> pos(g.number_of_vertices());
    size_t i = 0;
    for (size_t u = 0; u < g.number_of_vertices(); u++) {
      if (contains[u]) {
        g1.add_vertex(g.vertex_list[u]->vertex_data);
        pos[u] = i;
        i++;
      }
    }
    for (typename graph<V>::edge* e : g.edge_list) {
      size_t u = e->source_vertex->vertex_index;
      size_t v = e->target_vertex->vertex_index;
      if (contains[u] && contains[v]) {
        g1.add_edge(pos[u], pos[v], e->edge_weight);
      }
    }
    return g1;
  }

  template <class V>
  std::vector<digraph<V> > subgraphs(const digraph<V>& g, std::vector<size_t> component) {
    if (g.number_of_vertices() != component.size()) {
      throw std::invalid_argument("subgraph::vector size differs from number of vertices");
    }
    size_t max_component = 0;
    for (size_t i = 0; i < g.number_of_vertices(); i++) {
      if (component[i] > max_component) {
        max_component = component[i];
      }
    }
    std::vector<digraph<V> > g1(max_component + 1);
    std::vector<size_t> index(max_component + 1);
    std::vector<size_t> pos(g.number_of_vertices());
    for (size_t u = 0; u < g.number_of_vertices(); u++) {
      g1[component[u]].add_vertex(g.vertex_list[u]->vertex_data);
      pos[u] = index[component[u]];
      index[component[u]]++;
    }
    for (typename graph<V>::edge* e : g.edge_list) {
      size_t u = e->source_vertex->vertex_index;
      size_t v = e->target_vertex->vertex_index;
      if (component[u] == component[v]) {
        g1[component[u]].add_edge(pos[u], pos[v], e->edge_weight);
      }
    }
    return g1;
  }

  template <class V>
  digraph<V> transpose(const digraph<V>& g) {
    digraph<V> g1(g.number_of_vertices());
    for (size_t u = 0; u < g.number_of_vertices(); u++) {
      g1.vertex_list[u]->vertex_data = g.vertex_list[u]->vertex_data;
    }
    for (typename digraph<V>::edge* e : g.edge_list) {
      size_t u = e->source_vertex->vertex_index;
      size_t v = e->target_vertex->vertex_index;
      g1.add_edge(v, u, e->edge_weight);
    }
    return g1;
  }

}

#endif //DIGRAPH_HPP
