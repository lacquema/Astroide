# Create a python environment

A virtual environment is an isolated execution environment. It allows the packages used for a project to be isolated, so that there are different versions for each project.

To do this, use the following command: 

`python3 -m venv <environment_path>`

Note that the path to the python3 executable is in this environment: `<environment_path>/bin/python3`, and that the libraries are installed in the `<environment_path>/lib/` directory.


# Download the github directory

Now that the container is ready, you can download the gihub directory:

`git clone https://github.com/lacquema/Astroide`

# Compile and install packages

First, you need to open the newly installed github directory in your new environment:

`cd <environment_path>/Astroide/`

Now, you need to open the `<environment_path>/Astroide/Makefile` file and update the following parameters: COMPILF and PYTHON3. They correspond respectively to the paths to the executables of the fortran compiler installed on your computer and of python3 in your environment. The latter must therefore be `<environment_path>/bin/python3`.

And finally:

`make all`

This command compiles all the code and installs all the necessary packages.

# Launch the main code

`<environment_path>/bin/python3 <environment_path>/Astroide/Interface/Main.py`


