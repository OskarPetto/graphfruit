# graphfruit
**graphfruit** is a small graph library written in C++ that provides generic data structures and algorithms. It was created to learn basic C++ and to improve my knowledge in graph theory.

### Data structures:

* Graph
* Digraph
* Union-find
* Fibonacci heap

### Algorithms:

* Breadth-first search
* Cycle detection
* Dijkstra's shortest path
* Bellman-Ford shortest path
* Johnson's all shortest paths
* Kruskal's minimum spanning tree
* Prim's minimum spanning tree
* Khan's topological sorting

Output of *tester.cpp*:

```
-------Graph-------
2
4
-------Digraph-------
0
2
1
-------Breadth-first search-------
0
1
1
1
0
0
-------Shortest path-------
4
5
6
7
0
-------Shortest path (negative)-------
3
4
1
0
-------All shortest paths-------
0: 0 0 1 2 5 6 7 0 2
1: 0 1 1 2 5 2 7 1 2
2: 0 2 2 2 5 2 5 6 2
3: 0 2 3 3 3 2 5 6 2
4: 0 2 5 4 4 4 5 6 2
5: 0 2 5 2 5 5 5 6 2
6: 0 7 5 2 5 6 6 6 6
7: 0 7 5 2 5 6 7 7 7
8: 0 2 8 2 5 2 8 8 8
-------Is cyclic-------
0
-------Minimum spanning tree-------
(6, 7)
(8, 2)
(5, 6)
(2, 5)
(0, 1)
(2, 3)
(7, 0)
(3, 4)
-------Topological sorting-------
4
5
2
0
3
1
```
