/*
 * A class for directed graphs.
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
    explicit digraph(std::size_t n = 0);
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
    void add_edge(std::size_t source_vertex, std::size_t target_vertex, double edge_weight = 1.0);

    /*
     * Returns a vector of vectex pairs representing all edges.
     * Complexity: O(E)
     */
    std::vector<std::pair<std::size_t, std::size_t> > edges() const;

    /*
     * Returns the indegree of the vertex. Returns 0 if the vertex doesn't
     * exist.
     * Complexity: O(E)
     */
    std::size_t in_degree(std::size_t u) const;

    /*
     * Returns the number of directed edges in the graph.
     * Complexity: O(1)
     */
    std::size_t number_of_edges() const {return this->edge_list.size();}

    /*
     * Returns the outdegree of the vertex. Returns 0 if the vertex doesn't
     * exist.
     * Complexity: O(1)
     */
    std::size_t out_degree(std::size_t u) const;

    /*
     * Returns the predecessor vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(E)
     */
    std::vector<std::size_t> predecessors(std::size_t u) const;

    /*
     * Removes all edges between the source vertex and the target vertex from
     * the graph. Does nothing if there are no directed edges.
     * Complexity: O(E^2)
     */
    void remove_edges(std::size_t source_vertex, std::size_t target_vertex);

    /*
     * Returns the successor vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(V)
     */
    std::vector<std::size_t> successors(std::size_t u) const;

    /*
     * Uses Johnsons's algorithm to calculate the shortest paths between the
     * all pairs of vertices. Returns a 2D vector of vertices in the path in
     * reversed order. Returns an empty 2D vector if the graph contains a
     * negative cycle.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::vector<std::size_t> > johnson_all_shortest_paths(const digraph<V1>& g);

    /*
     * Uses Khan's algorithm to get a topological sorting of the graph. Returns
     * a vector of ordered vertices. If no such sorting is possible an empty
     * vector is returned.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<std::size_t> khan_topological_sort(digraph<V1>& g);

  };

  /*
   * digraph member implementations
   */
  template <class V>
  digraph<V>::digraph(std::size_t n)
    : base_graph<V>(n) {}

  template <class V>
  std::ostream& operator<<(std::ostream& out, const digraph<V>& g) {
    out << "Digraph: V=" << g.number_of_vertices();
    out << ", E=" << g.number_of_edges();
    for (typename digraph<V>::edge* e : g.edge_list) {
      out << std::endl;
      out << " ";
      out << "(" << e.source_vertex->vertex_index;
      out << ", " << e.target_vertex->vertex_index;
      out << ") " << e.edge_weight;
    }
    return out;
  }

  template <class V>
  void digraph<V>::add_edge(std::size_t source_vertex, std::size_t target_vertex, double edge_weight) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      std::size_t max = source_vertex > target_vertex ? source_vertex + 1 : target_vertex + 1;
      std::size_t i = this->number_of_vertices();
      this->vertex_list.resize(max);
      for (;i < max; i++) {
        this->vertex_list[i] = new vertex(i);
      }
    }
    edge* e = new edge(this->vertex_list[source_vertex], this->vertex_list[target_vertex], edge_weight);
    this->edge_list.push_back(e);
    this->vertex_list[source_vertex]->outgoing_edge_list.push_back(e); // seg fault
  }


  template <class V>
  std::vector<std::pair<std::size_t, std::size_t> > digraph<V>::edges() const {
    std::vector<std::pair<std::size_t, std::size_t> > v(number_of_edges());
    for (std::size_t i = 0; i < number_of_edges(); i++) {
      std::pair<std::size_t, std::size_t> a;
      a.first = this->edge_list[i]->source_vertex->vertex_index;
      a.second = this->edge_list[i]->target_vertex->vertex_index;
      v[i] = a;
    }
    return v;
  }

  template <class V>
  std::size_t digraph<V>::in_degree(std::size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    std::size_t d = 0;
    for (edge* e : this->edge_list) {
      if (e->target_vertex->vertex_index == u) {
        d++;
      }
    }
    return d;
  }

  template <class V>
  std::size_t digraph<V>::out_degree(std::size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    return this->vertex_list[u]->outgoing_edge_list.size();
  }

  template <class V>
  std::vector<std::size_t> digraph<V>::predecessors(std::size_t u) const {
    std::vector<std::size_t> v;
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
  void digraph<V>::remove_edges(std::size_t source_vertex, std::size_t target_vertex) {
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
  std::vector<std::size_t> digraph<V>::successors(std::size_t u) const {
    std::vector<std::size_t> v;
    if (!this->contains_vertex(u)) {
      return v;
    }
    for (edge* e : this->vertex_list[u].outgoing_edge_list) {
      v.push_back(e->target_vertex->vertex_index);
    }
    return v;
  }

  template <class V>
  std::vector<std::vector<std::size_t> > johnson_all_shortest_paths(const digraph<V>& g) {
    digraph<V> g1(g);
    std::size_t u = g1.number_of_vertices();
    g1.add_vertex();
    for (std::size_t v = 0; v < u; v++) {
      g1.add_edge(u, v, 0.0);
    }
    std::vector<double> bf_distance = g1.bellman_ford_distance(u);
    if (bf_distance.empty()) {
      std::vector<std::vector<std::size_t> > empty;
      return empty;
    }
    g1 = g;
    for (typename base_graph<V>::edge* e : g1.edge_list) {
      std::size_t u = e->source_vertex->vertex_index;
      std::size_t v = e->target_vertex->vertex_index;
      e->edge_weight += bf_distance[u] - bf_distance[v];
    }
    std::vector<std::vector<std::size_t> > previous(g1.number_of_vertices(), std::vector<std::size_t>(g1.number_of_vertices()));
    for (std::size_t i = 0; i < g1.number_of_vertices(); i++) {
      previous[i] = dijkstra_shortest_path(g1, i);
    }
    return previous;
  }

  template <class V>
  std::vector<std::size_t> khan_topological_sort(digraph<V>& g) {
    std::vector<std::size_t> sorted_vertices(g.number_of_vertices());
    std::vector<std::size_t> count(g.number_of_vertices());
    std::queue<std::size_t> q;
    for (typename digraph<V>::edge* e : g.edge_list) {
      count[e->target_vertex->vertex_index]++;
    }
    for (std::size_t i = 0; i < g.number_of_vertices(); i++) {
      if (count[i] == 0) {
        q.push(i);
      }
    }
    std::size_t i = 0;
    while (!q.empty()) {
      std::size_t u = q.front();
      q.pop();
      sorted_vertices[i] = u;
      i++;
      for (typename digraph<V>::edge* e : g.vertex_list[u]->outgoing_edge_list) {
        std::size_t v = e->target_vertex->vertex_index;
        count[v]--;
        if (count[v] == 0) {
          q.push(v);
        }
      }
    }
    if (i != g.number_of_vertices()) {
      std::vector<std::size_t> empty;
      return empty;
    }
    return sorted_vertices;;
  }

}

#endif //DIGRAPH_HPP
