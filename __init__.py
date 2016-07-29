from Tkinter import *
import scipy.cluster
from pymol import cmd
import numpy as np
import os
import sys
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import distutils.spawn
import MDAnalysis


path2 = os.path.abspath(os.path.dirname(__file__))
sys.path.append(path2)

####################################################################################################################################################
def __init__(self):

	self.menuBar.addmenuitem('Plugin', 'command', 'Apo-WaterDock', label = 'Apo-WaterDock', command = lambda s=self : option1())
	self.menuBar.addmenuitem('Plugin', 'command', 'Holo-Waterdock', label = 'Holo-waterdock', command = lambda s= self : option2())

###################################################################################################################################################
def option1():

	vinacomd = str(checkutilities())
	A = inputdata1()
	
	global proteinfile, ligandoption, ligandfile, centerx, centery, centerz

	if ligandoption == '2':

		UL = MDAnalysis.Universe(ligandfile)
		if ligandfile[-5:] == 'pdbqt':
			heavylig = UL.select_atoms('not type HD')
		else:
			heavylig = UL.select_atoms('not type H')

		com = np.mean(heavylig.positions, axis = 0)
		centerx = com[0]
		centery = com[1]
		centerz = com[2]


	waterfile()
	proteinpdbqtfile = proteinfile

	f1 = open('vinaconfig.txt','w')
	f1.write('receptor = ')
	f1.write(str(proteinpdbqtfile))
	f1.write('\nligand = water.pdbqt')
	f1.write('\nexhaustiveness = 20')
	f1.write('\nnum_modes = 20')
	f1.write('\ncenter_x = ')
	f1.write(str(centerx))
	f1.write('\ncenter_y = ')
	f1.write(str(centery))
	f1.write('\ncenter_z = ')
	f1.write(str(centerz))
	f1.write('\nsize_x = 15')
	f1.write('\nsize_y = 15')
	f1.write('\nsize_z = 15')
	f1.write('\nenergy_range = 100')
	f1.write('\nlog = outputlog.txt')
	f1.write('\nout = waterout.pdbqt')
	f1.close()

	for i in range(1,4):
		tempcomd = vinacomd + ' --config vinaconfig.txt'
		os.system(tempcomd)
		os.system("grep 'OW' waterout.pdbqt >> allwater.pdbqt")
		os.system("grep 'RESULT' waterout.pdbqt >> result.log")
		os.remove('outputlog.txt')
		os.remove('waterout.pdbqt')

	os.remove('vinaconfig.txt')
	coods = np.genfromtxt('allwater.pdbqt', usecols = (5,6,7), dtype = float)
	energy = np.genfromtxt('result.log', usecols = 3, dtype = float)

	selectcoods = np.compress(energy < -0.5, coods, axis = 0)
		
	fit1 = scipy.cluster.hierarchy.fclusterdata(selectcoods, 0.5, criterion = 'distance', metric = 'euclidean')
	fit1 = fit1.astype(int)
	numclust1 = np.max(fit1)

	clustcoods = np.zeros((numclust1,3), dtype = float)
	for i in xrange(1,numclust1+1):

		temp = np.compress(fit1 == i, selectcoods, axis = 0)
		tempavg = np.mean(temp, axis = 0)
		clustcoods[i-1,:] = tempavg

	fit2 = scipy.cluster.hierarchy.fclusterdata(clustcoods, 1.6, criterion = 'distance', metric = 'euclidean')
	fit2 = fit2.astype(int)
	numclust2 = np.max(fit2)

	finalcoods = np.zeros((numclust2,3), dtype = float)
	for i in xrange(1,numclust2+1):

		temp = np.compress(fit2 == i, clustcoods, axis = 0)
		tempavg = np.mean(temp, axis = 0)
		finalcoods[i-1,:] = tempavg

	if os.path.isfile('predictedwaters.pdb'):
		os.rename('predictedwaters.pdb','#predictedwaters.pdb')
	write_waterpdb('predictedwaters.pdb',finalcoods)

	os.remove('allwater.pdbqt')
	os.remove('result.log')
	os.remove('water.pdbqt')

	cmd.set("retain_order", 1)
	cmd.set("pdb_use_ter_records", 0)
	cmd.load(proteinfile, 'pro')
	cmd.show_as('cartoon', 'pro')
	cmd.color('blue','pro')
	cmd.load('predictedwaters.pdb','wats')
	cmd.color('red','wats')
	cmd.center('wats')
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
def option2():

	vinacomd = str(checkutilities())
	A = inputdata2()

	global proteinfile, ligandfile

	import addwater
	addwater.main(ligandfile)

	waterfile()
	#Writes the water.pdbqt file

	import dockcheck
	dockcheck.main(proteinfile, ligandfile, vinacomd)

	os.remove('waterdetails.txt')
	os.remove('placedwaters.pdb')
	os.remove('water.pdbqt')


	cmd.set("retain_order", 1)
	cmd.set("pdb_use_ter_records", 0)
	cmd.load(proteinfile, 'pro')
	cmd.show_as('cartoon', 'pro')
	cmd.color('blue','pro')
	cmd.load('predictedwaters.pdb','wats')
	cmd.color('red','wats')
	cmd.load(ligandfile,'lig')
	cmd.show_as('sticks','lig')
	cmd.center('lig')

#############################################################################################################################
def checkutilities():
	zz = distutils.spawn.find_executable('zina')
	homedir = str(os.environ.get('HOME'))
	file = homedir + '/pyvina.txt'
	if zz:
		if os.path.isfile(zz):
			vinacomd = zz

	elif os.path.isfile(file):
		fzz = open(file,'r')
		vinacomd = fzz.read()
		fzz.close()

	else:
		ZZ = vinapath()
		fzz = open(file,'r')
		vinacomd = fzz.read()
		fzz.close()

	return vinacomd
		
	
#############################################################################################################################
def waterfile():
	f1 = open('water.pdbqt', 'w')

	f1.write('REMARK  The pdbqt file for using water as a ligand\n')
	f1.write('ROOT\n')
	f1.write('ATOM      1  OW  HOH   231       0.950  11.375  16.494  1.00  0.00    -0.411 OA\n')
	f1.write('ATOM      2  HW1 HOH   231       1.766  11.375  17.071  1.00  0.00     0.205 HD\n')
	f1.write('ATOM      3  HW2 HOH   231       0.134  11.375  17.071  1.00  0.00     0.205 HD\n')
	f1.write('ENDROOT\n')
	f1.write('TORSDOF 0')

	f1.close()

#############################################################################################################################
def write_waterpdb(waterfilename, coordinates):
    
    numatom = coordinates.shape[0]
    xyz = open(waterfilename,'w')

    xyz.write("TITLE    Water Molecules Predicted by WaterDock\n")
    xyz.write("REMARK	Please Cite: Rapid and Accurate Prediction and Scoring of Water Molecules in Protein Binding Sites\n") 
    xyz.write("REMARK	DOI:10.1371/journal.pone.0032036\n")

    for i in xrange(0,numatom):
        header = 'HETATM'
        serial = i+1
        name = ' OW '
        resname = 'SOL'
        chainID = 'A'
        resSeq = i+1
        iCode = ' '
        occupancy = 1.0
        tempFactor = 0.0
        x, y, z = coordinates[i]
        xyz.write("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n" %(header,serial,name,iCode,resname,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor)) 

    xyz.close()
#############################################################################################################################
def unitvector(v):
	normal = np.linalg.norm(v)
	UV = v/normal
	return UV

#############################################################################################################################

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
class inputdata1:
	def __init__(self):

		global proteinfile

		self.window = Toplevel()
		self.window.resizable(0,0)
		self.window.title("Welcome to WaterDock 2.0")
		Label(self.window, text = "Predicting Waters in Apo-protein Structures").grid(row =1, sticky = E)

	
		self.ligopt = StringVar()
		self.profilename = StringVar()
		self.ligfilename = StringVar()
		self.xcom = StringVar()
		self.ycom = StringVar()
		self.zcom = StringVar()

		self.profilenameentry = StringVar()

		
		Label(self.window, text= "Import Protein from File (pdbqt)").grid(row=2, column=0, sticky=E)
		self.E1 = Entry(self.window, textvariable=self.profilenameentry)
		self.E1.grid(row=2, column=1)
		Button(self.window, text="Choose", command=self.profilechoose).grid(row=2, column=2, sticky = W)

	
		Label(self.window, text = "Identifying the binding site").grid(row=4,column=0, sticky = W)

		Radiobutton(self.window, text="Enter Center of Box (A)", variable = self.ligopt, value = '1').grid(row=5, column=0, sticky=E)
		Label(self.window, text='X').grid(row=5, column=1, sticky=E)
		Entry(self.window, textvariable=self.xcom).grid(row=5, column=2)
		Label(self.window, text='Y').grid(row=5, column=3, sticky=E)
		Entry(self.window, textvariable=self.ycom).grid(row=5, column=4)
		Label(self.window, text='Z').grid(row=5, column=5, sticky=E)
		Entry(self.window, textvariable=self.zcom).grid(row=5, column=6)

		Radiobutton(self.window, text="Co-ordinates from Ligand file (mol2/pdb/pdbqt)", variable = self.ligopt, value = '2').grid(row=6, column=0, sticky=E)
		self.E2 = Entry(self.window, textvariable=self.ligfilename)
		self.E2.grid(row=6, column=1)
		Button(self.window, text="Choose", command=self.ligfilechoose).grid(row=6, column=2, sticky = W)
	

		Button(self.window, text="Run", command=self.rungui).grid(row=13, column=3, sticky=E)
		Button(self.window, text="Cancel", command=self.byebye).grid(row=13, column=4, sticky=E)
		self.window.mainloop()

	def profilechoose(self):

		self.profilenamechoose = tkFileDialog.askopenfilename()

		self.E1.delete(0, END)
		self.E1.insert(0, self.profilenamechoose)


	def ligfilechoose(self):
		self.ligfilename = tkFileDialog.askopenfilename(filetypes=(("pdb files", "*.pdb"), ("pdbqt files", "*.pdbqt"), ("mol2 files", "*.mol2")))
		self.E2.insert(0, self.ligfilename)

	def byebye(self):
		self.window.destroy()

	def rungui(self):

		global proteinfile, ligandoption, ligandfile, centerx, centery, centerz
		if hasattr(self.ligopt, 'get'):
			ligandoption = self.ligopt.get()

		proteinfile = self.E1.get()

		if int(os.path.isfile(proteinfile)) == 0:
			tkMessageBox.showerror(title='Missing File', message='Protein File does not exist')
			sys.exit()
			
		if ligandoption == '2':
			if hasattr(self.E2, 'get'):
				ligandfile = self.E2.get()
				if int(os.path.isfile(ligandfile)) == 0:
					tkMessageBox.showerror(title='Missing File', message='Ligand File does not exist')
					sys.exit()

		elif ligandoption == '1':
			if hasattr(self.xcom,'get'):
				centerx = self.xcom.get()
				centerx = float(centerx)
			if hasattr(self.ycom,'get'):
				centery = self.ycom.get()
				centery = float(centery)
			if hasattr(self.zcom,'get'):
				centerz = self.zcom.get()
				centerz = float(centerz)

		if ligandoption == 2:

			UL = MDAnalysis.Universe(ligandfile)
			if ligandfile[-5:] == 'pdbqt':
				heavylig = UL.select_atoms('not type HD')
			else:
				heavylig = UL.select_atoms('not type H')

			com = np.mean(heavylig.positions, axis = 0)
			centerx = com[0]
			centery = com[1]
			centerz = com[2]

		self.window.quit()
########################################################################################################################
#############################################################################################################################
#############################################################################################################################
class inputdata2:
	def __init__(self):

		global proteinfile

		self.window = Toplevel()
		self.window.resizable(0,0)
		self.window.title("Welcome to WaterDock 2.0")
		Label(self.window, text = "Predicting Waters in Holo-protein Structures").grid(row =1, sticky = E)

	
		self.profilename = StringVar()
		self.ligfilename = StringVar()

		self.profilenameentry = StringVar()

		
		Label(self.window, text= "Import Protein from File (pdbqt)").grid(row=2, column=0, sticky=E)
		self.E1 = Entry(self.window, textvariable=self.profilenameentry)
		self.E1.grid(row=2, column=1)
		Button(self.window, text="Choose", command=self.profilechoose).grid(row=2, column=2, sticky = W)

	
		Label(self.window, text= "Import Ligand from File (pdb/mol2)").grid(row=3, column=0, sticky=E)
		self.E2 = Entry(self.window, textvariable=self.ligfilename)
		self.E2.grid(row=3, column=1)
		Button(self.window, text="Choose", command=self.ligfilechoose).grid(row=3, column=2, sticky = W)
	

		Button(self.window, text="Run", command=self.rungui).grid(row=4, column=3, sticky=E)
		Button(self.window, text="Cancel", command=self.byebye).grid(row=4, column=4, sticky=E)
		self.window.mainloop()

	def profilechoose(self):

		self.profilenamechoose = tkFileDialog.askopenfilename()

		self.E1.delete(0, END)
		self.E1.insert(0, self.profilenamechoose)


	def ligfilechoose(self):
		self.ligfilename = tkFileDialog.askopenfilename(filetypes=(("pdb files", "*.pdb"), ("pdbqt files", "*.pdbqt"), ("mol2 files", "*.mol2")))
		self.E2.insert(0, self.ligfilename)

	def byebye(self):
		self.window.destroy()

	def rungui(self):

		global proteinfile, ligandfile

		proteinfile = self.E1.get()


		if int(os.path.isfile(proteinfile)) == 0:
			tkMessageBox.showerror(title='Missing File', message='Protein File does not exist')
			sys.exit()
			
		if hasattr(self.E2, 'get'):
			ligandfile = self.E2.get()
			if int(os.path.isfile(ligandfile)) == 0:
				tkMessageBox.showerror(title='Missing File', message='Ligand File does not exist')
				sys.exit()

		self.window.quit()

########################################################################################################################
#############################################################################################################################
#############################################################################################################################

class vinapath:
	def __init__(self):

		self.window = Toplevel()
		self.window.resizable(0,0)
		self.window.title('Unable to find Vina executable')

		self.vinapathway = StringVar()

		Label(self.window, text = "Enter the path to the vina executable").grid(row = 1, column = 0, sticky = E)

		self.L1 = Entry(self.window, textvariable=self.vinapathway)
		self.L1.grid(row = 1, column = 1)

		Button(self.window, text="OK", command=self.okay).grid(row=2, column=1, sticky=E)
		self.window.mainloop()

	def okay(self):

		if hasattr(self.L1,'get'):
			self.vinapathway = self.L1.get()

		homedir = str(os.environ.get('HOME'))
		file = homedir + '/pyvina.txt'

		f1 = open(file,'w')
		f1.write(str(self.vinapathway))
		f1.close()

		self.window.quit()