# graphfruit
**graphfruit** is a small graph library written in C++. It was created to learn basic C++ and to improve my knowledge in graph theory. By using templates **graphfruit** provides generic data structures and algorithms.


Currently implemented algorithms:

* Breadth-first search
* Cycle detection
* Dijkstra's shortest path
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
