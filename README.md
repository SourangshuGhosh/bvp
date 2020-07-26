# bvp

Solution of Two-Point Boundary Value Problems
bvp, FENICS codes which use the finite element method to solve two point boundary value problems (BVP) over an interval in 1D.

Note that I have installed FENICS using Docker, and so to run this script I issue the commands:

```
cd $HOME/fenicsproject/bvp
fenicsproject run
python3 bvp_01.py
exit
```

Licensing:
The MIT license.

## Reference:

Hans Petter Langtangen, Anders Logg,
Solving PDEs in Python - The FEniCS Tutorial Volume 1.

## Source Code:

- bvp.sh runs all the tests
- bvp_01.py sets up and solves a two point boundary value problem.

- bvp_01.py the script.
- bvp_01.txt the output file.
- bvp_01_solution.png a plot of the solution.

**bvp_02.py repeats the calculation done by bvp_01.py, but now prints out a table of the solution, estimates the L2 norm and H1 seminorms of the error, and computes the average value of the solution in two ways.**

- bvp_02.py the script.
- bvp_02.txt the output file.

**bvp_03.py solves the same BVP 4 times, on a sequence of refined meshes.**

- bvp_03.py the script.
- bvp_03.txt the output file.
- bvp_03_04_mesh.xml a mesh.
- bvp_03_04_solution.png a plot of the corresponding solution.
- bvp_03_08_mesh.xml a mesh.
- bvp_03_08_solution.png a plot of the corresponding solution.
- bvp_03_16_mesh.xml a mesh.
- bvp_03_16_solution.png a plot of the corresponding solution.
- bvp_03_32_mesh.xml a mesh.
- bvp_03_32_solution.png a plot of the corresponding solution.

**bvp_04.py solves the same problem as bvp_01.py, but sets the Dirichlet boundary conditions explicitly.**

- bvp_04.py the script.
- bvp_04.txt the output file.
- bvp_04_solution.png a plot of the solution.

**bvp_05.py solves the same problem as bvp_04.py, except that on the right endpoint a Neumann condition is applied.**

- bvp_05.py the script.
- bvp_05.txt the output file.
- bvp_05_solution.png a plot of the solution.
