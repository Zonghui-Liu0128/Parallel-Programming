# Assignment 6: Molecular Dynamics(Part II)

## The Environment

- Processor: 2.70 GHz 2-Core Intel Core i5-5257U
- Logical cores/threads: 4
- OS: Windows 10
- Visual Studio Community 2022

## The Construct

- In `mian.c` file, I set the beginning of the parallel region at the outest loop. 
- In the outest loop, other functions should be run in single thread, so I set the `#pragma omp master`.
- In the `forces.c` file, I set `#pragma omp for` to use multiple threads to run the part.

## The Parameters

In `main.c` file:
- **Shared**: `f` and `x`, because they will be used outside.
- **firstprivate**: `side`, `sideh`, `rcoffs` and `npart`. Others will be initialised in the function.
- num_threads: `omp_set_num_threads(4);`

In `force.c` file:
- **reduction**: `epot` and `vir`

## Why Critical

- A critical section is a block of code which can be executed by only one thread at a time.
- In our situation, take the force as an example, it is the force between a partial and others. If we don't use critical, the shared variable will be overwritten and finally a wrong result.

## The Result

Sorry this part was not completed due to time limitation.

## Finding
Sorry this part was not completed due to time limitation.


















