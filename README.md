# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Installation and Compilation Guide

#### Installation
To get the files of this project, simply run this in your favorite terminal:

```git clone https://github.com/acse-2019/acse-4-sph-awe.git ```

#### Compilation
First, choose the files you want to implement. Usually user need to include .cpp and .h of certain method and SPH_Snippet,cpp file to run the project.

Choose **SPH_2D.cpp**and**SPH_2D.h** for forward Euler.

Choose **SPH_2D_pc.cpp**and**SPH_2D.h** for predictor_corrector.

Choose **SPH_2D_dynamic.cpp**and**SPH_2D_old_dynamic.h** for forward Euler with dynamic time stepping.

Choose **SPH_2D_pc_dynamic.cpp**and**SPH_2D_old_dynamic.h** for predictor_corrector with dynamic time stepping.

Notice that to implement another method you need to uncomment the corresponding functions in the header file and file_writter.cpp file.

For Visual Studio, please create a project under the home directory and only includes the files you need into the Header Files and Source Files folder. Then use the built-in compiler to compile these files.

For linux environment, simply run this on the command line:

```g++ fileYouNeedToUse.cpp fileYouNeedToUse.h -o SPH ```

It will generate an executable file in the root.

### User instructions

Our project contains the implementations of several methods for SPH.

#### For Windows:

After compling the code, user can run the executable file on the command line by inputting

``` SPH.exe ```

#### For Mac os and linuxs environment:

``` ./SPH ```

After running the the executable file, a simulation time output interval will be asked to enter. This interval will allow output files to be generated after a certain number of timesteps. So user will not get a huge number of files.

The output files are in the folder named **data** in the home directory by default. User can change this path in the *write_file* function of **file_writer.cpp** file.

#### For generating data of .vtp format:

The default file for generating data are **file_writer.cpp** and **file_writer.h**, whose output format of data is .vtp. And we provide a handful Jupyter notebook in the home directory to handle this data. By using this notebook user can get the data array of parameters like pressure, velocity etc, and the animation of the whole simulation process.

#### For generating data of .txt format:

We also provide the **data_output.cpp** and **data_output.h** file for generating the data of .txt format into the **post_process** folder(in the home directory), so user who is familiar with pandas can manipulate the data more easily. The Jupyter notebook named **txt_output.ipynb** in the home directory is an example of using .txt files to get the dataframe of all parameters.

### Documentation

We use **Doxygen** to generate the documentation file. Please see the index.html in the folder **html**.

Also we write a report for this project. Please see the **AWE_documentation**.

### Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
make runtests
```
under linux environment.
