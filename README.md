# WaterDock_pymol


[![PLOS ONE](https://img.shields.io/badge/PLOS%20ONE-10.1371%2Fjournal.pone.0172743-blue)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0172743)
[![PLOS ONE](https://img.shields.io/badge/PLOS%20ONE-10.1371%2Fjournal.pone.0032036-blue)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032036)

This is a PyMol Plugin for the original WaterDock and the new and improved WaterDock 2.0 specific for holo structures.

Written by Akshay Sridhar with updates from Patrick McCubbin and Phil Biggin.



Installation (linux and Mac but see below for Mac)
==================================================

There are 3 Files: __init__.py, addwater.py and dockcheck.py

Download all three and place in a folder and follow steps below. 

1) Open PyMol

2) Plugin --> Plugin Manager --> Install New Plugin --> Install From Local File

3) Choose the __init__.py file when prompted

4) DONE!!



Next time PyMol is restarted, Under Plugin, 2 Options will be available

1) Apo Waterdock -- The originial WaterDock

2) Holo Waterdock -- The new pipeline specific to holo structures



Dependencies

1) PyMol ofcourse 

2) MDAnalysis -- (version > 0.13)

3) numpy and Scipy (usually done as part of your MDanalysis install)

4) Autodock Vina 



Things to Note

1) The script runs AutoDock Vina through system calls of the command 'vina'. The system checks if the command 'vina' is callable and alerts if it is not. Ensure 'vina' in the system refers to AutoDock Vina. 

2) For the Holo-Dock, Please ensure the ligand has Hydrogens added to the structure. And inspect the correct protonation states of relevant atoms. 

3) The output file of both WaterDock options is written as 'predictedwaters.pdb'. If 'predictedwaters.pdb' already exists in the folder, it is renamed to 'predictedwaters1.pdb' and the new file is written as 'predictedwaters.pdb'.


For MacPymol Users
==================


M1 Apple Silicon Machines
=========================

Newer versions of MDAnalysis now support Apple Silicon (ie. M1 and M2 chips).  However, you will have to most likely compile vina from source (see below).  If you have vina working and installed and you already compiled pymol on your mac you can add then necessary packages thus (I write the packages as explicit separate lines here):

```
>>> import pip
>>> pip.main(['install','numpy'])
>>> pip.main(['install','scipy'])
>>> pip.main(['install','MDAnalysis'])
```

then follow the above installation WaterDock_Pymol

M1 Apple Silicon Machines - Vina
================================

At the time of writing (Oct 2023) there was pip installable version of vina for Apple Silicon.  Therefore, we have to compile ourselves

The instructions are here:-

 https://autodock-vina.readthedocs.io/en/latest/installation.html

But there are two things to note:

1.  You will need to use homebrew to install certain libraries.
2.  The Makefile will need changing to something like:-

```
BASE=/opt/homebrew
BOOST_VERSION=
BOOST_INCLUDE = $(BASE)/include
C_PLATFORM=-pthread
CXXFLAGS= DBOOST_TIMER_ENABLE_DEPRECATED
GPP=/usr/bin/clang++
C_OPTIONS= -O3 -DNDEBUG -std=c++11 -DBOOST_TIMER_ENABLE_DEPRECATED 
BOOST_LIB_VERSION=
include ../../makefile_common
```
Good luck!


Mac with intel chips
====================

Best installation route is via conda with python 3.5 --> 3.8 (we had trouble with python 3.9).

you can either compile pymol from source or the simplest route is pure conda with the schrondinger incentive verison of pymol as in the following example:-

```
conda create --name pymol-env python=3.7

conda activate pymol-env

conda config --add channels conda-forge

conda install  -c schrodinger pymol-bundle

conda install mdanalysis

conda install -c bioconda autodock-vina

conda install -c anaconda pyqt
```

Then fire up Pymol and install the plugin as per the instructions above.

Much Older Systems (legacy information only - probably only for systems older than 5 years)
============================================

1) Ensure that you have the correct version of python active

>> sudo port select --list python

This should return something like:-

>> Available versions for python:
>> none (active)
>> python26-apple
>> python27
>> python27-apple 

You need to make the python27-apple active:-

>> sudo port select python python27-apple

You can check it is by doing the list option again, but it should be fine, then fire up MacPymol and do the install and the plugin should work (assuming you have vina in your path as well)..
 

2) After Installation, copy the two additional files addwater.py and dockcheck.py into the folder you instructed pymol to install the plugin. 

Windows users
=============

Please use the separate windows-specific refactored branch.

 Help
 ====
 
 Please don't hesitate to report issues or email philip.biggin@bioch.ox.ac.uk
 
