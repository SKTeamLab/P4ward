.. Protacs Pipeline documentation master file, created by
   sphinx-quickstart on Thu Apr 20 13:59:38 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

P4ward: Predictive Protacs Python Pipeline - Workflow Automation for Research and Design.
===================================

P4ward is a fully automated, customizable and open-source program for Protacs ternary complex modelling. From two known binary complexes and a list of Protac 2D conformations, it can predict ternary complexes for each protac listed.

P4ward is written in Python3 and integrates permissively licensed structural biology tools, such as Megadock, RxDock, and libraries such as RDkit and BioPython to build a workflow with no accessibility barriers. As such, we acknowledge the countless contributors to open-source code which enables research by all and for all.

----------------
Getting Started
----------------

P4ward can be easily obtained with docker and conda. Please read :doc:`obtain` for more information.
For a gentle introduction to running P4ward, check the :doc:`tutorial`.

-----

.. toctree::
   :maxdepth: 1
   :caption: Getting Started
   
   obtain
   tutorial

.. toctree::
   :maxdepth: 1
   :caption: Documentation
   
   configure
   config_reference
   command_args
   modules_docs
