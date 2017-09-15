# graphfruit
**graphfruit** is a small graph library written in C++ that provides generic data structures and algorithms. It was created to learn basic C++ and to improve my knowledge in graph theory.

### Data structures:

* Graph
* Digraph
* Union-find
* Fibonacci heap

### Algorithms:

* Bellman-Ford shortest path
* Breadth-first search
* Cycle detection
* Dijkstra's shortest path
* Johnson's all shortest paths
* K-core decomposition
* Khan's topological sorting
* Kruskal's minimum spanning tree
* Longest path in DAGs
* Prim's minimum spanning tree

Output of *tester.cpp*:

```
-------Graph-------
2
4
-------Digraph-------
2
1
7
-------Breadth-first search-------
0
0
1
1
1
1
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
0: -1 0 1 2 5 6 7 0 2
1: 1 -1 1 2 5 2 7 1 2
2: 1 2 -1 2 5 2 5 6 2
3: 1 2 3 -1 3 2 5 6 2
4: 7 2 5 4 -1 4 5 6 2
5: 7 2 5 2 5 -1 5 6 2
6: 7 7 5 2 5 6 -1 6 6
7: 7 7 5 2 5 6 7 -1 7
8: 1 2 8 2 5 2 8 8 -1
-------Is cyclic and is DAG-------
0
1
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
0
1
2
3
4
5
-------Longest path in DAG-------
5
3
2
1
-------K-core-------
Graph: V=5, E=9
[0] - 3 - 2 - 1
[1] - 0 - 2 - 3 - 4
[2] - 0 - 3 - 1 - 4
[3] - 0 - 2 - 1 - 4
[4] - 1 - 2 - 3
```
