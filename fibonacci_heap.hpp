/*
 * A fibonacci heap implementation
 * @author Oskar Petto
 * @version 14.09.2017
 *  bug fixes
 * @version 13.09.2017
 *  first version
 */
#ifndef FIBONACCI_HEAP_HPP
#define FIBONACCI_HEAP_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

namespace graphfruit {

  /*
   * Forward declaration of fibonacci_heap class
   */
  template <class T, class Compare = std::less<T>>
  class fibonacci_heap;

  template <class T>
  class fibonacci_node {
    template <class T1, class Compare>
    friend class fibonacci_heap;
  public:

    /*
     * Constructors
     */
    explicit fibonacci_node(T key);
    fibonacci_node(const fibonacci_node& other) = default;
    fibonacci_node(fibonacci_node &&) noexcept = default;
    fibonacci_node& operator=(const fibonacci_node& other) = default;
    fibonacci_node& operator=(fibonacci_node&& other) noexcept = default;
    ~fibonacci_node() = default;

    T key;
    bool marked;
    std::size_t degree;
    fibonacci_node<T>* parent;
    fibonacci_node<T>* child;
    fibonacci_node<T>* left;
    fibonacci_node<T>* right;

  };

  template <class T, class Compare>
  class fibonacci_heap {
  public:

    /*
     * Constructors
     */
    explicit fibonacci_heap();
    fibonacci_heap(const fibonacci_heap& other);
    fibonacci_heap(fibonacci_heap &&) noexcept;
    const fibonacci_heap& operator=(const fibonacci_heap& other);
    const fibonacci_heap& operator=(fibonacci_heap&& other) noexcept;
    ~fibonacci_heap();

    /*
     * Makes the fibonacci heap printable to an output stream. The << operator
     * must be implemented by the type T.
     */
    template <class T1, class Compare1>
    friend std::ostream& operator<<(std::ostream& out, const fibonacci_heap<T1, Compare1>& h);

    /*
     * Removes all nodes in the heap and frees all used heap storage
     * Complexity: O(n)
     */
    void clear();

    /*
     * Decreases the key of node "node" to "key". If "key" is greater than the
     * key of "node", nothing happens.
     * Complexity: O(1) (amortized)
     */
    void decrease_key(fibonacci_node<T>* node, T key);

    /*
     * Returns true if and only if there are no nodes in the heap.
     * Complexity: O(1)
     */
    bool empty() const {return n == 0;}

    /*
     * Removes the min element of the heap.
     * Complexity: O(log(n)) (amortized)
     */
    void pop();

    /*
     * Inserts a new element in the heap and returns a pointer to the node.
     * Complexity: O(1)
     */
    fibonacci_node<T>* push(T key);

    /*
     * Returns the number of nodes in the heap.
     * Complexity: O(1)
     */
    std::size_t size() const {return n;}

    /*
     * Returns the top element of the heap.
     * Complexity: O(1)
     */
    T top() const;

  private:

    void clear(fibonacci_node<T>* node, std::size_t num);

    void consolidate();

    void cut(fibonacci_node<T>* node);

    void deep_copy(const fibonacci_heap& other);

    void deep_copy(fibonacci_node<T>* node);

    fibonacci_node<T>* link(fibonacci_node<T>* node1, fibonacci_node<T>* node2);

    void merge(fibonacci_node<T>* node1, fibonacci_node<T>* node2);

    std::ostream& print_level(std::ostream& out, fibonacci_node<T>* node, std::size_t level, bool after_first) const;

    void remove_min();

    void unlink(fibonacci_node<T>* node);

    Compare comp;
    std::size_t n;
    std::size_t trees;
    fibonacci_node<T>* min;

  };


  /*
   * fibonacci_node member implementations
   */
  template <class T>
  fibonacci_node<T>::fibonacci_node(T key)
    : key(key)
    , marked(false)
    , degree(0)
    , parent(nullptr)
    , child(nullptr)
    , left(this)
    , right(this) {}

  /*
   * fibonacci_heap member implementations
   */
  template <class T, class Compare>
  fibonacci_heap<T, Compare>::fibonacci_heap()
    : n(0)
    , trees(0)
    , min(nullptr)
    , comp(Compare()) {
    //empty
  }

  template <class T, class Compare>
  fibonacci_heap<T, Compare>::fibonacci_heap(const fibonacci_heap<T, Compare>& other) {
    deep_copy(other);
  }

  template <class T, class Compare>
  fibonacci_heap<T, Compare>::fibonacci_heap(fibonacci_heap<T, Compare>&& other) noexcept {
    deep_copy(other);
    other.clear();
  }

  template <class T, class Compare>
  const fibonacci_heap<T, Compare>& fibonacci_heap<T, Compare>::operator=(const fibonacci_heap<T, Compare>& other) {
    if (this == &other) {
      return *this;
    }
    clear();
    deep_copy(other);
    return *this;
  }

  template <class T, class Compare>
  const fibonacci_heap<T, Compare>& fibonacci_heap<T, Compare>::operator=(fibonacci_heap<T, Compare>&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    clear();
    deep_copy(other);
    other.clear();
    return *this;
  }

  template <class T, class Compare>
  fibonacci_heap<T, Compare>::~fibonacci_heap() {
    clear();
  }

  template <class T, class Compare>
  std::ostream& operator<<(std::ostream& out, const fibonacci_heap<T, Compare>& h) {
    if (h.empty()) {
      return out;
    }
    return h.print_level(out, h.min, 0, false);
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::clear() {
    clear(min, trees);
    n = 0;
    trees = 0;
    min = nullptr;
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::clear(fibonacci_node<T>* node, std::size_t num) {
    if (!node) {
      return;
    }
    for (std::size_t i = 0; i < num; i++) {
      fibonacci_node<T>* tmp = node->right;
      clear(node->child, node->degree);
      delete node;
      node = tmp;
    }
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::consolidate() {
    int max_degree = static_cast<int>(2.0 * log(static_cast<double>(n)));
    std::vector<fibonacci_node<T>*> degree_list(max_degree + 1);
    std::vector<fibonacci_node<T>*> root_list(trees);
    fibonacci_node<T>* root_node = min;
    for (std::size_t i = 0; i < trees; i++) {
      root_list[i] = root_node;
      root_node = root_node->right;
    }
    for (fibonacci_node<T>* node : root_list) {
      node->parent = nullptr;
      fibonacci_node<T>* tmp = node;
      while (degree_list[tmp->degree]) {
        tmp = link(tmp, degree_list[tmp->degree]);
        degree_list[tmp->degree-1] = nullptr;
      }
      degree_list[tmp->degree] = tmp;
    }
    for (fibonacci_node<T>* node : degree_list) {
      if (node && comp(node->key, min->key)) {
        min = node;
      }
    }
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::cut(fibonacci_node<T>* node) {
    if (!node->parent) {
      return;
    }
    if (node == node->right) {
      node->parent->child = nullptr;
    } else {
      node->parent->child = node->right;
    }
    unlink(node);
    merge(min, node);
    node->parent->degree--;
    node->parent = nullptr;
    node->marked = false;
    trees++;
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::decrease_key(fibonacci_node<T>* node, T key) {
    if (!node || comp(node->key, key)) {
      return;
    }
    node->key = key;
    if (node->parent && comp(node->key, node->parent->key)) {
      fibonacci_node<T>* parent = node->parent;
      cut(node);
      while (parent && parent->marked) {
        fibonacci_node<T>* tmp = parent->parent;
        cut(parent);
        parent = tmp;
      }
      if (parent) {
        parent->marked = true;
      }
    }
    if (comp(node->key, min->key)) {
      min = node;
    }
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::deep_copy(const fibonacci_heap& other) {
    if (other.empty()) {
      return;
    }
    deep_copy(other.min);
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::deep_copy(fibonacci_node<T>* node) {
    if (!node) {
      return;
    }
    fibonacci_node<T>* i = node;
    do {
      deep_copy(i->child);
      push(i->key);
      i = i->right;
    } while (i != node);
  }

  template <class T, class Compare>
  fibonacci_node<T>* fibonacci_heap<T, Compare>::link(fibonacci_node<T>* node1, fibonacci_node<T>* node2) {
    if (comp(node1->key, node2->key)) {
      fibonacci_node<T>* tmp = node1;
      node1 = node2;
      node2 = tmp;
    }
    unlink(node1);
    if (node2->child) {
      merge(node1, node2->child);
    } else {
      node2->child = node1;
    }
    node2->degree++;
    node1->parent = node2;
    node1->marked = false;
    trees--;
    return node2;
  }


  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::merge(fibonacci_node<T>* node1, fibonacci_node<T>* node2) {
    if (!node1 || !node2) {
      return;
    }
    node1->right->left = node2->left;
    node2->left->right = node1->right;
    node1->right = node2;
    node2->left = node1;
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::pop() {
    if (empty()) {
      return;
    }
    remove_min();
    if (!empty()) {
      consolidate();
    }
  }

  template <class T, class Compare>
  std::ostream& fibonacci_heap<T, Compare>::print_level(std::ostream& out,
    fibonacci_node<T>* node, std::size_t level, bool after_first) const {
    if (!node) {
      return out;
    }
    fibonacci_node<T>* i = node;
    do {
      if (after_first) {
        out << std::endl;
      }
      for (std::size_t j = 0; j < level; j++) {
        out << ' ';
      }
      out << i->key << ": " << i->degree;
      after_first = true;
      print_level(out, i->child, level + 1, after_first);
      i = i->right;
    } while (i != node);
    return out;
  }

  template <class T, class Compare>
  fibonacci_node<T>* fibonacci_heap<T, Compare>::push(T key) {
    fibonacci_node<T>* node = new fibonacci_node<T>(key);
    if (empty()) {
      min = node;
    } else {
      merge(node, min);
      if (comp(node->key, min->key)) {
        min = node;
      }
    }
    n++;
    trees++;
    return node;
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::remove_min() {
    fibonacci_node<T>* to_remove = min;
    if (to_remove == to_remove->right) {
      min = to_remove->child;
    } else {
      min = to_remove->right;
      if (to_remove->child) {
        merge(min, to_remove->child);
      }
    }
    unlink(to_remove);
    trees += to_remove->degree - 1;
    n--;
    delete to_remove;
  }

  template <class T, class Compare>
  T fibonacci_heap<T, Compare>::top() const {
    return min->key;
  }

  template <class T, class Compare>
  void fibonacci_heap<T, Compare>::unlink(fibonacci_node<T>* node) {
    if (!node) {
      return;
    }
    node->left->right = node->right;
    node->right->left = node->left;
    node->left = node;
    node->right = node;
  }

}

#endif //FIBONACCI_HEAP_HPP
