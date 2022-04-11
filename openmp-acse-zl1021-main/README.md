[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-f059dc9a6f8d3a56e377f745f24479a46679e63a5d9fe6f495e02850cd0d8118.svg)](https://classroom.github.com/online_ide?assignment_repo_id=6856350&assignment_repo_type=AssignmentRepo)
# Shared Memory Programming with OpenMP

## Course assignments

*Submission deadline: 5pm, 4th February 2022*

### Assignment 1: Hello World

This is a simple exercise to introduce you to the compilation and execution of
OpenMP programs. The example code can be found in `HelloWorld/`. First, compile
the code, making sure you use the appropriate flag to enable OpenMP. Before
running it, set the environment variable OMP_NUM_THREADS to a number `n`
between 1 and 4 with the command:

```
export OMP_NUM_THREADS=n
```

When run, the code enters a parallel region at the #pragma omp parallel
directive. At this point, `n` threads are spawned, and each thread executes the
print command separately.

The `omp_get_thread_num()` library routine returns a number (between `0` and
`n-1`) which identifies each thread. Incorporate a call to
`omp_get_num_threads()` into the code and print its value inside and outside
the parallel region.

### Assignment 2: Area of the Mandelbrot Set

This exercise demonstrates some of the issues that need to be considered when
adapting serial code to a parallel version.

*The Mandelbrot Set* is the set of complex numbers `c` for which the iteration:

```
z = z**2 + c
```

does not diverge from the initial condition `z=c`. To approximate whether a
point `c` lies in the set, a finite number of iterations are performed, and if
the condition `|z| > 2` is satisfied, then the point is considered to be
outside the Mandelbrot Set. What we are interested in is calculating the area
of the Mandelbrot Set. There is no known theoretical value for this, and
estimates are based on a procedure similar to that used here.

The method we shall use generates a grid of points in a box of the complex
plane containing the upper of the (symmetric) Mandelbrot Set. Then each point
is iterated using the equation above a finite number of times (say 2000). If
within that number of iterations the threshold condition `|z| > 2` is
satisfied, then that point is considered to be outside of the Mandelbrot Set.
Then counting the number of points within the Set and those outside will lead
to an estimate of the area of the Set.

Parallelise the serial code using the OpenMP directives and library routines
using the following method:

* Start a parallel region before the main loop, nest making sure that any private, shared or reduction variables within the region are correctly declared.
* Distribute the outermost loop across the threads available so that each thread has an equal number of the points. For this you will need to use some of the OpenMP library routines.

Once you have written the code try it out using 1, 2, 3 and 4 threads. Check
that the results are identical in each case, and compare the time taken for the
calculations using different number of threads. Is your solution well load
balanced? Try different ways of mapping iterations to threads.

Add a `README.md` file to the  Mandelbrot directory to describe your findings
in 200 words or less.
 
### Assignment 3: Mandelbrot revisited

Restart from the sequential code `Mandelbrot2/`. This time parallelise the
outer loop using a `parallel for` directive. Do not forget to use default(none)
and declare the shared, private and reduction variables. Add a schedule clause
and experiment with the different schedule kinds.

Add a `README.md` file to the  `Mandelbrot2` directory to describe your
findings in 200 words or less.

### Assignment 4: Nested parallelism

First, copy your code `Mandelbrot2` into a new folder `Mandelbrot3` (remember
to add this file to your repo). Now parallelise both outer loops in the
Mandelbrot example using nested parallel regions. To enable nested parallelism
use

```
export OMP_NESTED=true
```

and set the number of threads at each level with, for example,

```
export OMP_NUM_THREADS=2,4
```

You will need to think carefully about the data scoping (i.e.
shared/private/reduction). Experiment with different thread numbers – is the
performance ever better than with one level of parallelism? You can also try
parallelising both loops using a parallel loop directive on the outer loop and
a collapse clause.

Add a `README.md` file to the  `Mandelbrot3` directory to describe your
findings in 200 words or less.

### Assignment 5: Molecular Dynamics

The aim of this assignment is to explore the use of OpenMP critical constructs
to parallelise a molecular dynamics code.

The code can be found in `MolDyn/`. The code is a molecular dynamics (MD)
simulation of argon atoms in a box with periodic boundary conditions. The atoms
are initially arranged as a face-centred cubic (fcc) lattice and then allowed
to melt. The interaction of the particles is calculated using a Lennard-Jones
potential. The main loop of the program is in the file `main.c`.

Once the lattice has been generated and the forces and velocities initialised,
the main loop begins. The following steps are undertaken in each iteration of
this loop:

* The particles are moved based on their velocities, and the velocities are partially updated (call to `domove`).
* The forces on the particles in their new positions are calculated and the virial and potential energies accumulated (call to `forces`).
* The forces are scaled, the velocity update is completed and the kinetic energy calculated (call to `mkekin`).
* The average particle velocity is calculated and the temperature scaled (call to `velavg`).
* The full potential and virial energies are calculated and printed out (call to `prnout`).

#### Parallelisation method

The parallelisation of this code is a little less straightforward. There are
several dependencies within the program which will require use of the critical
construct as well as the reduction clause. The instructions for parallelising
the code are as follows:

Edit the function in `forces.c`. Add an `omp parallel` directive to the
outermost loop in this subroutine, identifying any private or reduction
variables. Hint: There are two reduction variables.  Identify the variable
within the loop which must be protected from possible race conditions and use
the `omp critical` directive to ensure this is the case.

Once this is done the code should be ready to run in parallel. Compare the
output using 2, 3 and 4 threads with the serial output to check that it is
working.

Try adding the schedule clause with the kind `static,n` to the `for` directive
for different values of `n`. Does this have any effect on performance?

Add a `README.md` file to the  `MolDyn` directory to describe your findings in
200 words or less.

### Assignment 6: Molecular Dynamics Part II

Copy your code from `MolDyn` in the previous assignment to a new folder
`MolDyn` (remember to add these files to your repo).

Following the previous assignment, we will update the molecular dynamics code
to take advantage of orphaning and examine the performance issues of using
critical regions.

To reduce the overhead of starting and stopping the threads, you can change the
`parallel for` directive to a `for` directive and start the parallel region
outside the program's main loop. In `main.c` enclose the main loop in a
parallel region. Except for the forces routine, all the other work in this loop
should be executed by one thread only. Ensure that this is the case using the
`single` or `master` construct.

Recall that any reduction variables updated in a parallel region but outside of
the `for` construct should be updated by one thread only. As before, check your
code is still working correctly.

Is there any difference in performance compared with the code without
orphaning?

Try using array reductions, or atomic directives, or lock routines (with one
lock variable per particle) instead of critical regions. How did this impact
the parallel performance of the code?

Add a `README.md` file to the  `MolDyn2` directory to describe your findings in
200 words or less.
