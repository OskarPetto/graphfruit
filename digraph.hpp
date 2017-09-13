/*
 * A class for directed graphs.
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
     * Adds a directed edge between source_vertex and target_vertex to
     * the graph and adds the vertices if they don't exist to the graph.
     * Complexity: O(1)
     */
    void add_edge(std::size_t source_vertex, std::size_t target_vertex, double edge_weight = 1.0);

    /*
     * Returns true if and only if there exists a directed edge between
     * source_vertex and target_vertex. Returns false if either vertex is not
     * in this graph.
     * Complexity: O(V)
     */
    bool contains_edge(std::size_t source_vertex, std::size_t target_vertex) const;

    /*
     * Returns the indegree of the vertex. Returns 0 if the vertex
     * does not exist.
     * Complexity: O(E)
     */
    std::size_t in_degree(std::size_t u) const;

    /*
     * Returns the outdegree of the vertex. Returns 0 if the vertex
     * does not exist.
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
     * Removes all edges between the source_vertex and the target_vertex from
     * the graph. Does nothing if there are no directed edges.
     * Complexity: O(E * E)
     */
    void remove_edges(std::size_t source_vertex, std::size_t target_vertex);

    /*
     * Returns the successor vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(V)
     */
    std::vector<std::size_t> successors(std::size_t u) const;

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
    for (typename base_graph<V>::edge* e : g.edge_list) {
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
  bool digraph<V>::contains_edge(std::size_t source_vertex, std::size_t target_vertex) const {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      return false;
    }
    for (edge* e : this->vertex_list[source_vertex]->outgoing_edge_list) {
      if (e->target_vertex == this->vertex_list[target_vertex]) {
        return true;
      }
    }
    return false;
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
