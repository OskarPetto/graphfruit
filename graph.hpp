/*
 * A class for undirected graphs.
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
    explicit graph(std::size_t n = 0);
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
     * Returns the degree of the vertex. Returns 0 if the vertex does not exist.
     * Complexity: O(1)
     */
    std::size_t degree(std::size_t u) const;

    /*
     * Returns the neighbour vertices of vertex. Return an empty vector if
     * the vertex doesn't exist.
     * Complexity: O(E)
     */
    std::vector<std::size_t> neighbours(std::size_t u) const;

    /*
     * Removes all edges between the source_vertex and the target_vertex from
     * the graph. Does nothing if there are no directed edges.
     * Complexity: O(E * E)
     */
    void remove_edges(std::size_t source_vertex, std::size_t target_vertex);

    /*
     * Uses Kruskal's algorithm to find the edges of a minimum spanning tree of
     * the graph. Returns the edges a vector of vertex pairs. If the graph is
     * not connnected the result is not defined.
     * Complexity: O(E * log(V))
     */
    template <class V1>
    friend std::vector<std::pair<std::size_t, std::size_t>> kruskal_minimum_spanning_tree(graph<V1>& g);

    /*
     * Uses Prims's algorithm to find the edges of a minimum spanning tree of
     * the graph. Returns the edges a vector of vertex pairs. If the graph is
     * not connnected the result is not defined.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::pair<std::size_t, std::size_t>> prim_minimum_spanning_tree(graph<V1>& g);

  protected:

    // for sort in kruskal_minimum_spanning_tree
    struct edge_comp {
      bool operator() (edge* a, edge*b) {
        return a->edge_weight <= b->edge_weight;
      }
    };

  };

  /*
   * graph member implementations
   */
  template <class V>
  graph<V>::graph(std::size_t n)
    : base_graph<V>(n) {}

  template <class V>
  std::ostream& operator<<(std::ostream& out, const graph<V>& g) {
    out << "Graph: V=" << g.number_of_vertices();
    out << ", E=" << g.number_of_edges();
    for (typename base_graph<V>::edge* e : g.edge_list) {
      out << std::endl;
      out << " ";
      out << "{" << e.source_vertex->vertex_index;
      out << ", " << e.target_vertex->vertex_index;
      out << "} " << e.edge_weight;
    }
    return out;
  }

  template <class V>
  void graph<V>::add_edge(std::size_t source_vertex, std::size_t target_vertex, double edge_weight) {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      std::size_t max = source_vertex > target_vertex ? source_vertex + 1: target_vertex + 1;
      std::size_t i = this->number_of_vertices();
      this->vertex_list.resize(max);
      for (;i < max; i++) {
        this->vertex_list[i] = new vertex(i);
      }
    }
    edge* e1 = new edge(this->vertex_list[source_vertex], this->vertex_list[target_vertex], edge_weight);
    edge* e2 = new edge(this->vertex_list[target_vertex], this->vertex_list[source_vertex], edge_weight);
    this->edge_list.push_back(e1);
    this->vertex_list[source_vertex]->outgoing_edge_list.push_back(e1);
    this->vertex_list[target_vertex]->outgoing_edge_list.push_back(e2);
  }

  template <class V>
  bool graph<V>::contains_edge(std::size_t source_vertex, std::size_t target_vertex) const {
    if (!this->contains_vertex(source_vertex) || !this->contains_vertex(target_vertex)) {
      return false;
    }
    if (this->vertex_list[source_vertex]->outgoing_edge_list.size() <
        this->vertex_list[target_vertex]->outgoing_edge_list.size()) {
      for (edge* e : this->vertex_list[source_vertex]->outgoing_edge_list) {
        if (e->target_vertex == this->vertex_list[target_vertex]) {
          return true;
        }
      }
    } else {
      for (edge* e : this->vertex_list[target_vertex]->outgoing_edge_list) {
        if (e->target_vertex == this->vertex_list[source_vertex]) {
          return true;
        }
      }
    }
    return false;
  }

  template <class V>
  std::size_t graph<V>::degree(std::size_t u) const {
    if (!this->contains_vertex(u)) {
      return 0;
    }
    return this->vertex_list[u]->outgoing_edge_list.size();
  }

  template <class V>
  std::vector<std::size_t> graph<V>::neighbours(std::size_t u) const {
    std::vector<std::size_t> v;
    if (!this->contains_vertex(u)) {
      return v;
    }
    for (edge* e : this->vertex_list[u]->outgoing_edge_list) {
      v.push_back(e->target_vertex->vertex_index);
    }
    return v;
  }

  template <class V>
  void graph<V>::remove_edges(std::size_t source_vertex, std::size_t target_vertex) {
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
  std::vector<std::pair<std::size_t, std::size_t>> kruskal_minimum_spanning_tree(graph<V>& g) {

    std::vector<std::pair<std::size_t, std::size_t>> result(g.number_of_vertices() - 1);
    std::vector<typename graph<V>::edge*> sorted_edges(g.edge_list.begin(), g.edge_list.end());
    std::sort(sorted_edges.begin(), sorted_edges.end(), typename graph<V>::edge_comp());
    union_find uf(g.number_of_vertices());
    std::size_t i = 0;
    for (typename graph<V>::edge* e : sorted_edges) {
      if (i == g.number_of_vertices() - 1) {
        return result;
      }
      std::size_t u_root = uf.set_find(e->source_vertex->vertex_index);
      std::size_t v_root = uf.set_find(e->target_vertex->vertex_index);
      if (u_root != v_root) {
        result[i] = std::make_pair(e->source_vertex->vertex_index, e->target_vertex->vertex_index);
        i++;
        uf.set_union(u_root, v_root);
      }
    }
    return result;
  }

  template <class V>
  std::vector<std::pair<std::size_t, std::size_t>> prim_minimum_spanning_tree(graph<V>& g) {

    std::vector<std::pair<std::size_t, std::size_t>> result(g.number_of_vertices() - 1);
    std::vector<double> distance(g.number_of_vertices());
    std::vector<fibonacci_node<std::pair<std::size_t, double>>*> fib_nodes(g.number_of_vertices());
    fibonacci_heap<std::pair<std::size_t, double>, typename base_graph<V>::fib_comp> min_heap;

    for (std::size_t i = 0; i < g.number_of_vertices(); i++) {
      distance[i] = std::numeric_limits<double>::max();
      std::pair<std::size_t, double> a(i, std::numeric_limits<double>::max());
      fib_nodes[i] = min_heap.push(a);
    }
    distance[0] = 0.0;
    min_heap.decrease_key(fib_nodes[0], std::make_pair(0, 0.0));
    while (!min_heap.empty()) {
      std::pair<std::size_t, double> u = min_heap.top();
      fib_nodes[u.first] = nullptr;
      min_heap.pop();
      for (typename base_graph<V>::edge* e : g.vertex_list[u.first]->outgoing_edge_list) {
        std::size_t v = e->target_vertex->vertex_index;
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
