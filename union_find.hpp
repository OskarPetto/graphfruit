/*
 * A union find implementation using weighted quick-union and path compression.
 * @author Oskar Petto
 * @version 13.09.2017
 *  first version
 */
#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <cstddef>
#include <vector>

namespace graphfruit {

  class union_find {
  public:

    /*
     * Constructors
     */
    explicit union_find(std::size_t n);
    union_find(const union_find& other) = default;
    union_find(union_find &&) noexcept = default;
    union_find& operator=(const union_find& other) = default;
    union_find& operator=(union_find&& other) = default;
    ~union_find() = default;

    /*
     * Removes all elements from the structure.
     * Complexity: O(1)
     */
    void clear();

    /*
     * Uses one-pass path compression to find the root of the set containing x.
     * Complexity: O(log(n))
     */
    std::size_t set_find(std::size_t x);

    /*
     * Uses union by rank to merge the sets of x and y.
     * Complexity:
     */
    void set_union(std::size_t x, std::size_t y);

  private:

    struct elem {
      std::size_t parent;
      unsigned int rank;
    };

    std::vector<elem> id;

  };

}

#endif //UNION_FIND_HPP
