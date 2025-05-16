# Table of contents:

- [Installation on personal computeur](#install_on_computer)

- [Installation on computing server](#install_on_server)

- [Simulation](#simulation)
    - [Choice of integrator](#integrator)
    - [Generation of sub-simulations](#generation)
    - [Launch](#launch)
    - [Continuation](#continuation)
    - [Extraction](#extraction)

- [Analyse](#analyse)

<br><br>

<div id='install_on_computer'>  

# Installation on personal computer

## Prerequisites

You need python3 and a fortran compiler, as gfortran.

## Create a container

It is recommended that you create a virtual environment. It is an isolated execution environment. It allows the packages used for a project to be isolated, so that there are different versions for each project.

To do this, use the following command: 

`python3 -m venv <environment_path>`

You can now activate the environment: 

`source <environment_path>/bin/activate`

The `python3` command will now refer to the python3 installed in your environment `<environment_path>/bin/python3` and not the one on the system. And the libraries are now installed in the `<environment_path>/lib/` directory.

If you do not choose this option, simply create an usual directory:

`mkdir <environment_path>`


## Download the github directory  

Now that the container is ready, you can download the gihub directory inside:

`cd <environment_path>`

`git clone https://github.com/lacquema/Astroide`


## Install python packages

All the necessary python packages are listed in the `<environment_path>/Astroide/requirements.txt` file. You can install them all:

`python3 -m <environment_path>/bin/pip install -r <environment_path>/Astroide/requirements.txt`


## Compile fortran code

If you are using the fortran compiler `gfortran`, you can directly compile all the fortran files:

`cd <environment_path>/Astroide`

`make compile`

Else, open the `<environment_path>/Astroide/Makefile` file and update the following parameters: COMPILF. This corresponds to the paths or simply the command for the executables of the fortran compiler installed on your computer. Then, run the commands above. 

<br><br>

<div id='install_on_server'>  

# Installation on computing server

You can import the code from github in the same way as on your personal computer. The difference lies in the installation of the software dependencies. 

On computing servers, the use of modern package managers such as Guix is highly recommended and often mandatory. These tools enable isolated and reproducible management of software dependencies, thus avoiding unintentional changes to the server environment. In addition, Guix avoids unnecessary duplication of packages by storing shared dependencies only once in the server's `/gnu/store`. This considerably reduces the use of storage space. 

Instead of installing packages directly, you configure them to reference the Guix-managed library. All required packages are explicitly listed in the `<environment_path>/Astroide/manifest.scm` file. You can reference them all:

`guix package -m <environment_path>/Astroide/manifest.scm`

Ensure to consult the Guix documentation provided by your server administrators, as specific configurations or permissions might apply.

Note that software dependencies may includes the python packages we need, as well as python3 itself and the fortran compiler `gfortran-toolchain`.

At this point, all that remains is to compile the fortran files in the same way as on a personal computer. 


<br><br>

<div id='simulation'>  

# Simulation

<div id='integrator'>  

## Choice of integrator

This code currently offers a choice of two integrators: RMVS3 and HJS, which differ in the way they calculate system dynamics.

RMVS3 uses heliocentric coordinates for each body. It is well suited to planetary systems in which the planets have negligible mass compared to the central star.

HJS, on the other hand, is based on Jacobi coordinates for each body in the system. It is more expensive in terms of computation time, but offers better accuracy than RMVS3. This integrator is particularly well suited to hierarchical multiple star systems.

<div id='generation'>  

## Generation of sub-simulations

This code divides the work into sub-simulations to optimize computing time. Indeed, it's not the integration of planet dynamics that requires the most resources, but that of test particles, or debris, which is often far more numerous. To remedy this, the code decouples the calculation of particle dynamics, spreading them over several independent sub-simulations. This enables artificial parallelization: each sub-simulation can be run separately, in parallel, which significantly reduces total computation time.

Once the integrator has been selected, you need to refer to the subsimulation generation file, available in the following folder :

`<environment_path>/Astroide/Generator`

Two scripts are available to generate sub-simulations, depending on the integrator chosen: `gen_multi_rmvs3.sh` for RMVS3 and `gen_multi_hjs.sh` for HJS.

Both scripts require input parameters specific to the simulation context. The special feature of HJS is its use of Jacobi coordinates, which means that the hierarchy of orbits must be specified: for each body, you must indicate whether the others are located within its orbit, or whether they are indifferent to its dynamics at this scale. There are therefore two generation shell files, one for RMVS3 and another for HJS.

<div id='launch'>  

## Launch

<div id='continuation'>  

## Continuation

<div id='extraction'>  

## Extraction



<br><br>

<div id='analyse'>  

# Analyse
    `python3 <environment_path>/Astroide/Interface/Main.py`