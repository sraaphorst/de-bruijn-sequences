# de Bruijn Cycles

Simple Python script to generate [de Bruijn cycles](https://en.wikipedia.org/wiki/De_Bruijn_sequence)
for parameters `n` and `k`, where:
* `n` is the size of the each word; and
* `k` is the alphabet size.

The aim is to have, as a final result, a cycle that covers all `k^n` possible words in the alphabet.

Here is an example of a simple cycle where `n = 3` and `k = 2`:

```text
$ ./de-bruijn.py 2 3
[0, 0, 0, 1, 0, 1, 1, 1]
```

We can see that the output cycles to cover all `2^3 = 8` words:

<table>
  <thead>
    <tr>
      <th colspan="2">Codewords</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>000</td>
      <td>011</td>
    </tr>
    <tr>
      <td>001</td>
      <td>111</td>
    </tr>
    <tr>
      <td>010</td>
      <td>110</td>
    </tr>
    <tr>
      <td>101</td>
      <td>100</td>
    </tr>
  </tbody>
</table>



This works by:
* Creating the digraph of all `k^(n-1)` vertices.
* Creating edges from all vertices:
  * `v = (t_0, t_1, ..., t_{n-2})` to
  * `w = (t_1, ..., t_{n-2}, s)`.
* Ensuring that the digraph has the necessary properties:
  * The out degree of each vertex is the same as its in degree.
  * The digraph forms a strongly connected component (there is always a path between any two vertices). This is done with [Kosaraju's algorithm](https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm).
* Then an Eulerian circuit is found using [Hierholzer's Algorithm](https://en.wikipedia.org/wiki/Eulerian_path).
* The final cycle is normalized (any trailing zero is moved to the beginning) and checked for coverage.