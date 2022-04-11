# The Structure of the HPC Files

The structure tree is showed in the tree.txt file.

- **CMakeLists.txt**: the configuration file used by `Cmake` to build `Makefile`.
- **Makefile**: is used in by the `make` tool to control the compilation and linking of engineering projects.
- **object_files**: the running results from the HPC PBS system, just like the results we get on our local laptop after the instruction of `mpiexec`.
- **scripts**: the scripts used for the HPC PBS system.
- **build**: all the files(size = 100 * 100)
- **build3**: all the files(size = 300 * 300)
- **build5**: all the files(size = 500 * 500)
- **build10**: all the files(size = 1000 * 1000)

# Instructions
In the root path,

- Upload my source code file and the CMakeLists.txt
- Create a new folder: `mkdir build`
- Enter into the folder: `cd build`
- Use CMake tool to generate the Makefile: `cmake ..`
- Build the code using make tool: `make`
- Submit the job using PBS system: `qsub script100-1.pbs`
- Check the status of all of the jobs: `qstat -a`






















