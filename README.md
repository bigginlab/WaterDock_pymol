# WaterDock-2.0

Please Cite:
Rapid and Accurate Prediction and Scoring of Water Molecules in Protein Binding Sites

DOI : http://dx.doi.org/10.1371/journal.pone.0032036

This is a PyMol Plugin for the original WaterDock and the new and improved WaterDock 2.0 specific for holo structures.

INSTALLATION

There are 3 File. __init__.py, addwater.py and dockcheck.py

Download all three and place in a folder and follow steps below. 

1) Open PyMol

2) Plugin --> Plugin Manager --> Install New Plugin --> Install From Local File

3) Choose the __init.py file when prompted

4) DONE!!



Next time PyMol is restarted, Under Plugin, 2 Options will be available

1) Apo Waterdock -- The originial WaterDock

2) Holo Waterdock -- The new pipeline specific to holo structures



Python Dependencies

1) PyMol ofcourse

2) MDAnalysis -- (version > 0.13)

3) numpy -- (any version compatible with the MDAnalysis build)

4) scipy



Things to Note

1) The script runs AutoDock Vina through system calls of the command 'vina'. The system checks if the command 'vina' is callable and alerts if it is not. Ensure 'vina' in the system refers to AutoDock Vina. 

2) For the Holo-Dock, Please ensure the ligand has Hydrogens added to the structure. And inspect the correct protonation states of relevant atoms. 

3) The output file of both WaterDock options is written as 'predictedwaters.pdb'. If 'predictedwaters.pdb' already exists in the folder, it is renamed to 'predictedwaters1.pdb' and the new file is written as 'predictedwaters.pdb'. 
