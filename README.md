# ACOTSP - Local Search

A study of Ant Colony Optimization with Local Search on the Travelling Salesman Problem. This project was done for the COMP 4520 Honours Project course at the University of Manitoba under the supervision of Dr. Parimala Thulasiraman.

Program execution is defined as follows:

```(bash)
Usage: ./acotsp -i in_file
ACO parameters:
  -a: weight of pheromone on decision
  -b: weight of heuristic data on decision
  -g: max number of generations
  -l: use local search (2-opt) number of iterations
  -N: number of ants
  -Q: amount of pheromone to be deposited
  -q: degree of random choice
  -r: rate of pheromone evaporation
  -s: seed for random choice
  -t: initial pheromone levels
  -n: number of ACO runs
  -h: help
Input/Output:
  -i: specify tsp file to use as input/source file
  -o: specify file to output results to
```
