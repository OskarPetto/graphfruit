#include "union_find.hpp"

namespace graphfruit {

  union_find::union_find(std::size_t n) {
    id.resize(n);
    for (std::size_t i = 0; i < n; i++) {
      id[i].parent = i;
      id[i].rank = 1;
    }
  }

  void union_find::clear() {
    id.clear();
  }

  std::size_t union_find::set_find(std::size_t x) {
    while (id[x].parent != x) {
      id[x].parent = id[id[x].parent].parent; // one-pass
      x = id[x].parent;
    }
    return x;
  }

  void union_find::set_union(std::size_t x, std::size_t y) {
    std::size_t x_root = set_find(x);
    std::size_t y_root = set_find(y);
    if (x_root == y_root) {
      return;
    }
    if (id[x_root].rank < id[y_root].rank) {
      id[x_root].parent = y_root;
    } else if (id[x_root].rank > id[y_root].rank) {
      id[y_root].parent = x_root;
    } else {
      id[y_root].parent = x_root;
      id[x_root].rank++;
    }
  }

}
