# Coursework 

Deadline: 25 February 2022 at 5pm

Initial files - Coursework description and serial code


## Description of Folders/Files

- `tree.txt`: the structure tree for the project.
- `src`: the source code for the project, including the cpp file and the CMakeLists.txt file.
- `post-processing`: the files for the post-processing part, more details can be find in the README.md file in that folder.
- `HPCFiles`: the object files and script files for the HPC system, more details can be find in the README.md file in that folder. 

ps: I apologize that due to time limitation, I upload all the files from HPC files to the repo. Actually. I have submitted many jobs and only some of them have been used for figures.

- `fig`: the pictures used in the report.
- `report`: the report of the project


## Parameters

the parameters in the `my_code.cpp` file in the src folder:

- To set the number of points in the whole domain: `int imax = 100, jmax = 100;`
- To set the size of the whole domain: `double x_max = 10.0, y_max = 10.0;`
- To set the iteration time: `double t_max = 30.0;`
- Use the Neumann boundary conditions: `bool is_Neumann = true;` and `bool is_periodic = false;`
- Use the Dirichlet boundary conditions: `bool is_Neumann = false;` and `bool is_periodic = false;`
- Use the Periodic boundary conditions: `bool is_Neumann = false;` and `bool is_periodic = true;`


## Instructions
### Local
- **Build**: mpicxx -o my_code my_code.cpp
- **Run**: mpiexec -n 6 my_code

### HPC
- mkdir build
- cd build
- cmake ..
- make
- qsub my_scripts.pbs
- qstat -a

More details could be seen in the `HPCFiles` folder

## Reference
[1] https://www.cnblogs.com/zhouzhe-blog/p/9614360.html

[2] https://blog.csdn.net/sinat_29957455/article/details/85145488

[3] Lecture2 worksheet solutions

[4] Lecture4 worksheet solutions
