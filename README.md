# Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure

[![Documentation Status](https://readthedocs.org/projects/dg-wave-c/badge/?version=latest)](https://dg-wave-c.readthedocs.io/en/latest/?badge=latest)

<!--ts-->
   * [Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure](#dynamic-load-balancing-for-a-hp-adaptive-discontinuous-galerkin-wave-equation-solver-via-space-filling-curve-and-advanced-data-structure)
      * [Introduction](#introduction)
      * [Setup](#setup)
         * [What You Need](#what-you-need)
         * [Compile and Execute](#compile-and-execute)
      * [Documentation](#documentation)
      * [Source Code Documentation](#source-code-documentation)
      * [Approximation of Wave Equation](#approximation-of-wave-equation)
      * [AMR Refinement Types: hp-adaptivity](#amr-refinement-types-hp-adaptivity)
         * [h-refinement](#h-refinement)
         * [p-refinement](#p-refinement)
      * [AMR Data Structure: Hash Table](#amr-data-structure-hash-table)
      * [Dynamic Load Balancing](#dynamic-load-balancing)
         * [Traditional Way: Graph-based Repartitioning Algorithm](#traditional-way-graph-based-repartitioning-algorithm)
         * [My Proposal: Space-Filling Curves (SFCs) Based Repartitioning Algorithm](#my-proposal-space-filling-curves-sfcs-based-repartitioning-algorithm)
      * [Program Proformace](#program-proformace)
         * [Speedup](#speedup)

<!-- Added by: shiqi, at: Wed Dec  2 20:33:38 EST 2020 -->

<!--te-->

## Introduction
We combine a high-order method -- the **discontinuous Galerkin [spectral element method](https://en.wikipedia.org/wiki/Spectral_element_method) (DG-SEM)**, 
with parallel [**adaptive mesh refinement and coarsening (AMR)**](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) techniques and apply it to a **two-dimensional [wave equation](https://en.wikipedia.org/wiki/Wave_equation) solver**.

Advanced data structures and dynamic load balancing are applied to the solver to achieve efficient data management and high-level parallelism. 

## Setup
### What You Need
* C++
* [CMake](https://cmake.org/) (at least version 3.9)
* [GCC](https://gcc.gnu.org/) 7.5.0 (GNU Compiler Collection)
* [OpenMPI](https://www.open-mpi.org/) 4.0.2

### Compile and Execute
```
mkdir build && cd build
cmake ..
make 
mpirun -np 4 main
```
You also can use parallel `make`, *e.g.*, `make -j 16`. 16 threads will be used to compile the code. 

`mpirun -np 4 main` executes the program with 4 processors, you could change the processor number as you want.

## Documentation
A detailed documentation an be found at [here](https://dg-wave-c.readthedocs.io/en/latest/).

## Source Code Documentation
Source code explanation:
[Source code documentation]( https://shiqihe000.github.io/DG_wave_c/doxygen/html/index.html)

## Approximation of Wave Equation
The basic model of wave propagation is the wave equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" title="\frac{\partial ^{2}p}{\partial t^{2}}-c^{2}(p_{xx}+p_{yy})=0" /></a>

The variable `p` represents the acoustic pressure and `c` is the sound speed. 

## AMR Refinement Types: hp-adaptivity
Two types of refinements are implemented in this work: h-refinement and p-refinement. 
### h-refinement
<p align="center">
  <img src="./imgs/h_refinement.png" width="100" height = "40" >
</p>
Subdivide an element into children elements. 

### p-refinement
<p align="center">
  <img src="./imgs/p_refinement.png" width="100" height = "40" >
</p>
Raise polynomial orders inside the targeted element. 

## AMR Data Structure: Hash Table
In stead of using the conventional [tree data structure](https://en.wikipedia.org/wiki/Tree_(data_structure)), 
[hash table data structure](https://en.wikipedia.org/wiki/Hash_table) is employeed. 

## Dynamic Load Balancing
The adaptivity of the program introduces **load imbalance** among the processors, causing the program to slow down. 
In order to achieve high-level parallelism, dynamic load balancing strategy is apllied to balance the workload. 

The workload redistribution problem is an optimization problem with **two main goals**: 
the workload should be distributed evenly among the processors with small memory overhead and the
interfacing boundaries between partitions should be as small as possible. The optimization
problem has been proven to be **NP-hard**.

### Traditional Way: Graph-based Repartitioning Algorithm
pros:
* Well-studied.
* Multiple libraries.

cons:
* Reached scalability limits (requires "global" graph knowledge). 
* High memory consumption. 

### My Proposal: Space-Filling Curves (SFCs) Based Repartitioning Algorithm
pros:
* Simplifies a multi-dimensional partitioning problem into a one-dimensional one. 
* SFCs have good locality and can be generated fast. 
* Low memory usage. 
* Favours distributed systems. 

cons: 
* The partitioning boundaries are decided by the trajectry of the space-filling curve. 

[Hilber curve](https://en.wikipedia.org/wiki/Hilbert_curve) is chosen for this work. 
The figure below shows the first three levels of Hilbert curve. 

<p align="center">
  <img src="./imgs/Hilbert_curve2.png" width="400" height = "120" >
</p>


## Program Proformace
### Speedup

<p align="center">
  <img src="./imgs/speedup_factor_cedar.png" width="500" height = "300" >
</p>

The test is done on [Cedar](https://docs.computecanada.ca/wiki/Cedar) with processor number 
ranging from 32 to 2048. We gained a maximum speedup around 8. 
