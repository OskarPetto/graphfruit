/*
 * An abstract class for undirected and directed graphs.
 * @version 13.09.2017
 *  first version
 */

#ifndef BASE_GRAPH_HPP
#define BASE_GRAPH_HPP

#include <cstddef>
#include "fibonacci_heap.hpp"
#include <iostream>
#include <queue>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace graphfruit {

  /*
   * This is an abstract class which is inherited by graph and digraph.
   * V is the type of data that can be stored in a vertex.
   */
  template <class V = void*>
  class base_graph {
  public:

    /*
     * Constructors
     */
    explicit base_graph(std::size_t n = 0);
    base_graph(const base_graph<V>& other);
    base_graph(base_graph<V>&& other) noexcept;
    base_graph& operator=(const base_graph<V>& other);
    base_graph& operator=(base_graph<V>&& other) noexcept;
    ~base_graph();

    /*
     * Operators to compare two graphs.
     */
    bool operator==(const base_graph<V>& other) const;
    bool operator!=(const base_graph<V>& other) const;

    /*
     * Adds an edge between source vertex and target_vertex to the graph and
     * adds vertices to the graph if they don't exist.
     * Complexity: O(1) amortized
     */
    virtual void add_edge(std::size_t source_vertex, std::size_t target_vertex, double edge_weight = 1.0) = 0;

    /*
     * Creates a new vertex and adds it to the graph.
     * Complexity: O(1) amortized
     */
    void add_vertex(V vertex_data = V());

    /*
     * Removes all vertices and edges and reinitializes an empty graph.
     * Complexity: O(V + E)
     */
    void clear();

    /*
     * Returns true if and only if there exists an edge between source vertex
     * and target vertex. Returns false if either vertex is not in this graph.
     * Complexity: O(V)
     */
    virtual bool contains_edge(std::size_t source_vertex, std::size_t target_vertex) const = 0;

    /*
     * Returns true if and only if there exists a vertex with the given index.
     * Complexity: O(1)
     */
    bool contains_vertex(std::size_t u) const;

    /*
     * Copies all vertices and edges from the other graph.
     * Complexity: O(V + E)
     */
    void deep_copy(const base_graph<V>& other);

    /*
     * Returns the weight of the edge. Throws an illegal argument exception if
     * the either vertex or the edge doesn't exist.
     * Complexity: O(V)
     */
    double edge_weight(std::size_t source_vertex, std::size_t target_vertex) const;

    /*
     * Returns a vector of size_t pairs representing the source and target
     * vertices of all edges.
     * Complexity: O(E)
     */
    std::vector<std::pair<std::size_t, std::size_t>> edges() const;

    /*
     * Returns true if and only if there are no vertices in the graph.
     * Complexity: O(1)
     */
    bool empty() const {return vertex_list.size() == 0;}

    /*
     * Returns the number of edges in the graph.
     * Complexity: O(1)
     */
    std::size_t number_of_edges() const {return edge_list.size();}

    /*
     * Returns the number of vertices in the graph.
     * Complexity: O(1)
     */
    std::size_t number_of_vertices() const {return vertex_list.size();}

    /*
     * Removes all edges between the source vertex and the target vertex from
     * the graph. Does nothing if there are no edges.
     * Complexity: O(E * E)
     */
    virtual void remove_edges(std::size_t source_vertex, std::size_t target_vertex) = 0;

    /*
     * Removes a vertex from the graph. Removing a vertex also removes all
     * edges that have this vertex as source or as target. Does nothing if the
     * vertex does not exist.
     * Complexity: O(V * V + E * E)
     */
    void remove_vertex(std::size_t u);

    /*
     * Returns the data stored in the vertex. Throws an illegal argument
     * exception if the vertex doesn't exist.
     * Complexity: O(1)
     */
    V vertex_data(std::size_t u) const;

    /*
     * Performs a BFS on the graph starting from start vertex.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend std::vector<bool> breadth_first_search(const base_graph<V1>& g, std::size_t start_vertex);

    /*
     * Uses Dijkstra's algorithm to calculate the shortest paths between the
     * start vertex and all other vertices. Returns a vector of predecessors in
     * these shortest paths. Returns a vector with all zeros if start vertex is
     * not in the graph.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::size_t> dijkstra_shortest_path(const base_graph<V1>& g, std::size_t start_vertex);

    /*
     * Uses Dijkstra's algorithm to calculate the shortest paths between the
     * start vertex and the end vertex. Returns a vector of vertices in the
     * path in reversed order. Returns an empty vector if the start vertex or
     * the end vertex are not in the graph.
     * Complexity: O(E + V * log(V))
     */
    template <class V1>
    friend std::vector<std::size_t> dijkstra_shortest_path(const base_graph<V1>& g, std::size_t start_vertex, std::size_t end_vertex);

    /*
     * Uses a modified DFS to detect cycles in a graph. Returns true if and
     * only if the graph contains a cycle.
     * Complexity: O(V + E)
     */
    template <class V1>
    friend bool is_cyclic(const base_graph<V1>& g);

  protected:

    struct vertex;

    struct edge {

      vertex* source_vertex; //the tail vertex
      vertex* target_vertex; //the head vertex
      double edge_weight; // cost/distance of the edge

      explicit edge(vertex* source_vertex, vertex* target_vertex, double edge_weight = 1.0)
        : source_vertex(source_vertex)
        , target_vertex(target_vertex)
        , edge_weight(edge_weight) {}

    };

    struct vertex {

      std::size_t vertex_index; // index of vertex in graph
      V vertex_data; // store data in this vertex
      std::vector<edge*> outgoing_edge_list; // edges leaving this vertex

      explicit vertex(std::size_t vertex_index, V vertex_data = V())
        : vertex_index(vertex_index)
        , vertex_data(vertex_data) {}

    };

    std::vector<vertex*> vertex_list;
    std::vector<edge*> edge_list;

    // for fibonacci_heap in dijkstra_shortest_path
    struct fib_comp {
      bool operator() (const std::pair<std::size_t, double>& a, const std::pair<std::size_t, double>& b) {
        return a.second < b.second;
      }
    };

    bool is_cyclic_util(std::vector<int>& color, std::size_t u) const;

  };

  /*
   * graph member implementations
   */
  template <class V>
  base_graph<V>::base_graph(std::size_t n) {
    vertex_list.resize(n);
    for (std::size_t i = 0; i < n; i++) {
      vertex_list[i] = new vertex(i);
    }
  }

  template <class V>
  base_graph<V>::base_graph(const base_graph<V>& other) {
    deep_copy(other);
  }

  template <class V>
  base_graph<V>::base_graph(base_graph<V>&& other) noexcept {
    deep_copy(other);
    other.clear();
  }

  template <class V>
  base_graph<V>& base_graph<V>::operator=(const base_graph<V>& other) {
    if (this == &other) {
      return *this;
    }
    clear();
    deep_copy(other);
    return *this;
  }

  template <class V>
  base_graph<V>& base_graph<V>::operator=(base_graph<V>&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    clear();
    deep_copy(other);
    other.clear();
    return *this;
  }

  template <class V>
  base_graph<V>::~base_graph() {
    clear();
  }

  template <class V>
  bool base_graph<V>::operator==(const base_graph<V>& other) const {
    //TODO: isomorphism
    return false;
  }

  template <class V>
  bool base_graph<V>::operator!=(const base_graph<V>& other) const {
    return !(*this == other);
  }

  template <class V>
  void base_graph<V>::add_vertex(V vertex_data) {
    vertex_list.push_back(new vertex(vertex_list.size(), vertex_data));
  }

  template <class V>
  void base_graph<V>::clear() {
    for (edge* e : edge_list) {
      delete e;
    }
    for (vertex* u: vertex_list) {
      delete u;
    }
    edge_list.clear();
    vertex_list.clear();
  }

  template <class V>
  bool base_graph<V>::contains_vertex(std::size_t u) const {
    return u < vertex_list.size();
  }

  template <class V>
  void base_graph<V>::deep_copy(const base_graph<V>& other) {
    std::size_t n = vertex_list.size();
    for (vertex* u : other.vertex_list) {
      add_vertex(u->vertex_data);
    }
    for (edge* e : other.edge_list) {
      vertex* u = e->source_vertex;
      vertex* v = e->target_vertex;
      add_edge(u->vertex_index + n, v->vertex_index + n, e->edge_weight);
    }
  }

  template <class V>
  double base_graph<V>::edge_weight(std::size_t source_vertex, std::size_t target_vertex) const {
    if (!contains_vertex(source_vertex) || !contains_vertex(target_vertex)) {
      throw std::invalid_argument("edge_weight::vertex does not exist");
    }
    for (edge* e : vertex_list[source_vertex]->outgoing_edge_list) {
      if (e->target_vertex == vertex_list[target_vertex]) {
        return e->edge_weight;
      }
    }
    throw std::invalid_argument("edge_weight::edge does not exist");
  }

  template <class V>
  std::vector<std::pair<std::size_t, std::size_t>> base_graph<V>::edges() const {
    std::vector<std::pair<std::size_t, std::size_t>> v(number_of_edges());
    for (std::size_t i = 0; i < number_of_edges(); i++) {
      std::pair<std::size_t, std::size_t> a;
      a.first = edge_list[i]->source_vertex->vertex_index;
      a.second = edge_list[i]->target_vertex->vertex_index;
      v[i] = a;
    }
    return v;
  }

  template <class V>
  void base_graph<V>::remove_vertex(std::size_t u) {
    if(!contains_vertex(u)) {
      return;
    }
    typename std::vector<edge*>::iterator i = edge_list.begin();
    while (i != edge_list.end()) {
      if ((*i)->source_vertex == vertex_list[u] || (*i)->target_vertex == vertex_list[u]) {
        (*i)->target_vertex->in_degree--;
        delete *i;
        i = edge_list.erase(i);
      } else {
        ++i;
      }
    }
    vertex_list[u]->outgoing_edge_list.clear();
    for (vertex* v : vertex_list) {
      typename std::vector<edge*>::iterator j = v->outgoing_edge_list.begin();
      while (j != v->outgoing_edge_list.end()) {
        if ((*j)->target_vertex == vertex_list[u]) {
          j = v->outgoing_edge_list.erase(j);
        } else {
          ++j;
        }
      }
    }
    typename std::vector<vertex*>::iterator k = vertex_list.begin();
    while(k != vertex_list.end()) {
      if ((*k)->vertex_index == u) {
        delete *k;
        k = vertex_list.erase(k);
      } else {
        if ((*k)->vertex_index > u) {
          (*k)->vertex_index--;
        }
        ++k;
      }
    }
  }

  template <class V>
  V base_graph<V>::vertex_data(std::size_t u) const {
    if (!contains_vertex(u)) {
      throw std::invalid_argument("vertex_data::vertex does not exist");
    }
    return vertex_list[u]->vertex_data;
  }

  template <class V>
  std::vector<bool> breadth_first_search(const base_graph<V>& g, std::size_t start_vertex) {
    std::vector<bool> processed(g.number_of_vertices());
    if (!g.contains_vertex(start_vertex)) {
      return processed;
    }
    processed[start_vertex] = true;
    std::queue<std::size_t> q;
    q.push(start_vertex);
    while (!q.empty()) {
      std::size_t u = q.front();
      q.pop();
      for (typename base_graph<V>::edge* e : g.vertex_list[u]->outgoing_edge_list) {
        std::size_t v = e->target_vertex->vertex_index;
        if (!processed[v]) {
          processed[v] = true;
          q.push(v);
        }
      }
    }
    return processed;
  }

  template <class V>
  std::vector<std::size_t> dijkstra_shortest_path(const base_graph<V>& g, std::size_t start_vertex) {
    std::vector<std::size_t> previous(g.number_of_vertices());
    if (!g.contains_vertex(start_vertex)) {
      return previous;
    }
    std::vector<double> distance(g.number_of_vertices());
    std::vector<fibonacci_node<std::pair<std::size_t, double>>*> fib_nodes(g.number_of_vertices());
    fibonacci_heap<std::pair<std::size_t, double>, typename base_graph<V>::fib_comp> min_heap;

    for (std::size_t i = 0; i < g.number_of_vertices(); i++) {
      distance[i] = std::numeric_limits<double>::max();
      std::pair<std::size_t, double> a(i, std::numeric_limits<double>::max());
      fib_nodes[i] = min_heap.push(a);
    }
    distance[start_vertex] = 0.0;
    min_heap.decrease_key(fib_nodes[start_vertex], std::make_pair(start_vertex, 0.0));
    previous[start_vertex] = start_vertex;

    while (!min_heap.empty()) {
      std::pair<std::size_t, double> u = min_heap.top();
      fib_nodes[u.first] = nullptr;
      min_heap.pop();
      for (typename base_graph<V>::edge* e : g.vertex_list[u.first]->outgoing_edge_list) {
        std::size_t v = e->target_vertex->vertex_index;
        if (fib_nodes[v] && distance[v] > distance[u.first] + e->edge_weight) {
          distance[v] = distance[u.first] + e->edge_weight;
          min_heap.decrease_key(fib_nodes[v], std::make_pair(v, distance[v]));
          previous[v] = u.first;
        }
      }
    }
    return previous;
  }

  template <class V>
  std::vector<std::size_t> dijkstra_shortest_path(const base_graph<V>& g, std::size_t start_vertex, std::size_t end_vertex) {
    std::vector<std::size_t> path;
    if (!g.contains_vertex(start_vertex) || !g.contains_vertex(end_vertex)) {
      return path;
    }
    std::vector<std::size_t> previous = dijkstra_shortest_path(g, start_vertex);
    std::size_t i = end_vertex;
    while (i != start_vertex) {
      path.push_back(i);
      i = previous[i];
    }
    path.push_back(start_vertex);
    return path;
  }

  template <class V>
  bool is_cyclic(const base_graph<V>& g) {
    if (g.empty()) {
      return false;
    }
    std::vector<int> color(g.number_of_vertices());
    for (typename base_graph<V>::vertex* u : g.vertex_list) {
      if (color[u->vertex_index] == 0) {
        if (g.is_cyclic_util(color, u->vertex_index)) {
          return true;
        }
      }
    }
    return false;
  }

  template <class V>
  bool base_graph<V>::is_cyclic_util(std::vector<int>& color, std::size_t u) const {
    color[u] = 1;
    for (edge* e : vertex_list[u]->outgoing_edge_list) {
      std::size_t v = e->target_vertex->vertex_index;
      if (color[v] == 1) {
        return true;
      }
      if (color[v] == 0 && is_cyclic_util(color, v)) {
        return true;
      }
    }
    color[u] = 2;
    return false;
  }

}

#endif // BASE_GRAPH_HPP
