- **animation_Dirichlet.gif**: the animation when using Dirichlet boundary conditions(set the flag `is_Neumann = false`).
- **animation_Neumann.gif**: the animation when using Neumann boundary conditions(set the flag `is_Neumann = true`).
- **Compare_performance.xlsx**: the running time in `*.o` file from HPC system. In this file, I also calculate the speedup ratio and the parallel efficiency and  then generate line graphs to compare performance.
- **post-processing.ipynb**: the python code for combining the status files from each processor and generating the animation.

ps: the code of generating animation is from the internet