Obtaining P4ward
================

The easiest way to obtain and run P4ward is through Docker, where it is possible to obtain and run P4ward through a single command. It is also possible to run the code using a conda environment as well as using Apptainer (former singularity) for running in high performance computational clusters (HPCs).

Docker
------

Running using a single command by obtaining it from Docker Hub:

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

Install Megadock for single node environment, CPU (`Instructions here <https://github.com/akiyamalab/MEGADOCK/blob/master/doc/BUILD.md#d-compile-for-cpu-node-only-thread-parallelization>`_).

Clone the P4ward github repo:

.. code:: bash

    git clone https://github.com/PaulaJLR/p4ward.git
    
Then build a conda environment to work P4ward:

.. code:: bash

    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda create -n protacs_pipeline python=3.11 --file ./p4ward/dockerfiles/conda_requirements.txt

Next, add P4ward to ``PYTHONPATH``:

.. code:: bash

    export PYTHONPATH=$PYTHONPATH:path/to/cloned/p4ward/

When running P4ward with this strategy, always remember to activate the conda environment, add the program to ``PYTHONPATH``, and then run:

.. code:: bash

    conda activate p4ward
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

After the ``.sif`` file has been generated, P4ward can be run by:

.. code:: bash

    apptainer run -B /[root_mount_path] /path/to/p4ward.sif --config config.ini

The root mount path will be whichever filesystem you are working on in the cluster, for example, it could be ``-B /scratch`` or ``-B /project``.