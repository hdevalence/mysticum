mysticum
========

software for displaying Pascal's mysticum

notes
-----

* The sage code is in `src/sage/mysticum.sage`
* Run unit tests with `sage mysticum-test.sage`

todo
----

*   more efficient algorithms for working with permutations
*   add tests for steiner nodes
*   add tests for cayley lines
*   refactor code to have a mysticum object that keeps all of the
    data of the points and lines in the configuration, rather than 
    recomputing all the data all the time
*   have a good way to represent "permutation up to inverse" so
    that we don't need to worry about choosing a representative
    correctly / consistently
 

