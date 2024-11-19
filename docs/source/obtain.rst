Obtaining P4ward
================

The easiest way to obtain and run P4ward is through Docker, where it is possible to obtain and run P4ward through a single command. It is also possible to run the code using a conda environment as well as using Apptainer (former singularity) for running in high performance computational clusters (HPCs).

Docker
------

Running using a single command by obtaining it from Docker Hub (this option will be available soon):

.. code:: bash

    sudo docker run -v .:/home/data paulajlr/p4ward:latest --config_file config.ini

This will obtain the container, if it hasn't been done already, and directly run the pipeline.

----------

If you have obtained the docker container through a tar file, such as ``p4ward.tar``, you can build the container with:

.. code:: bash

    sudo docker load -i path/to/p4ward.tar

Then the program can be run using the container with:

.. code:: bash

    sudo docker run -v .:/home/data p4ward --config_file config.ini

----------

You may also choose to build the container yourself. In this case, clone the P4ward git repository to your chosen working directory, navigate into the new ``p4ward`` folder and build:

.. code:: bash

    git clone https://github.com/PaulaJLR/p4ward.git
    cd p4ward
    sudo docker build -f dockerfiles/Dockerfile -t p4ward

Then run as previously described:

.. code:: bash

    sudo docker run -v .:/home/data p4ward --config_file config.ini


Conda
-----

It is also possible to easily obtain P4ward's dependencies and run the program locally without the need for containers. This step assumes you have conda or miniconda installed.

Install Megadock for single node environment, CPU (`Instructions here <https://github.com/akiyamalab/MEGADOCK/blob/master/doc/BUILD.md#d-compile-for-cpu-node-only-thread-parallelization>`_). Make sure to add megadock installation to your system PATH, with:

.. code:: bash

    export PATH=$PATH:path/to/megadock/installation/

Test if this step was successful by navigating to another directory and calling megadock with ``megadock -h``. If you want this PATH update to be permanent, you can paste it on your ``~/.bashrc`` file. Otherwise, remember to add megadock to PATH before running P4ward.

Clone the P4ward github repo:

.. code:: bash

    git clone https://github.com/SKTeamLab/P4ward.git
    
Then build a conda environment to work P4ward:

.. code:: bash

    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda create -n p4ward python=3.11 --file ./p4ward/dockerfiles/conda_requirements.txt

Next, add P4ward to ``PYTHONPATH``:

.. code:: bash

    export PYTHONPATH=$PYTHONPATH:path/to/cloned/p4ward/

.. tip::

   More about PYTHONPATH:
   When you run python from the command line, there is a way to tell your system where to find additional python packages, so that you can call them without having to worry about their installation path in your system. When you obtain the program, you get the main folder for the repository, which is called "P4ward". Inside of it, there are many resources. For example, there is a folder called "tutorial", where you can find all the files for this tutorial, another called "docs", where this documentation was written, and, importantly, another folder called P4ward, which contains the python package itself. This means that python should find this folder, and so you must add the root P4ward repo folder to PYTHONPATH. This way, python will look in the root p4ward repo folder and in it, it will find the "P4ward" package. Thus, if for example you ran ``git clone https://github.com/SKTeamLab/P4ward.git`` while you were in your Downloads folder, you should add this path to your PYTHONPATH: ``/home/USER/Downloads/p4ward``
    

When running P4ward with this strategy, always remember to add megadock to your PATH, activate the conda environment, add the program to ``PYTHONPATH``, and then run the program. Example:

.. code:: bash

    conda activate p4ward
    export PATH=$PATH:path/to/megadock/installation/
    export PYTHONPATH=$PYTHONPATH:path/to/cloned/p4ward/
    python -m p4ward --config_file config.ini


Apptainer
---------

Usually, HPC clusters do not support Conda or Docker. In this case, it is possible to convert a Docker container into an Apptainer container, which is usually supported by clusters. If you don't have the tar file of the Docker container, you can make one by running:

.. code:: bash

    sudo docker save -o p4ward.tar p4ward

Next, a Docker tarfile can be converted to an Apptainer file by running:

.. code:: bash

    apptainer build p4ward.sif docker-archive:p4ward.tar

At this point, make sure to do the conversion within a job script. It should not need more than 16GB of RAM and not much more than an hour to run, depending on the node specifications. This was tested with Apptainer version 1.3.4 and above.
After the ``.sif`` file has been generated, P4ward can be run by:

.. code:: bash

    apptainer run -B /[root_mount_path] /path/to/p4ward.sif --config config.ini

The root mount path will be whichever filesystem you are working on in the cluster, for example, it could be ``-B /scratch`` or ``-B /project``.