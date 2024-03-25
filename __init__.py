from pymol.Qt import QtWidgets as qtw
from pymol.Qt import QtCore as qtc
from pymol import cmd
from pymol import plugins
import numpy as np
import scipy.cluster
import os
import sys
import distutils.spawn
import MDAnalysis


path2 = os.path.abspath(os.path.dirname(__file__))
sys.path.append(path2)


def __init_plugin__(app=None):
    plugins.addmenuitemqt('Apo-WaterDock', option1)
    plugins.addmenuitemqt('Holo-waterdock', option2)


dialog = None  # global reference avoids garbage collection


def option1():

    global dialog
    vinacomd = str(checkutilities())
    dialog = inputdata1()

    if dialog.result() == 1:
        runapowaterdock(vinacomd, dialog.proteinfile, dialog.centerx, dialog.centery, dialog.centerz)


def runapowaterdock(vinacomd, proteinfile, centerx, centery, centerz):

    waterfile()

    f1 = open('vinaconfig.txt', 'w')
    f1.write('receptor = ')
    f1.write(str(proteinfile))
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
    f1.write('\nverbosity = 1')
    f1.write('\nout = waterout.pdbqt')
    f1.close()

    for i in range(1, 4):
        tempcomd = vinacomd + ' --config vinaconfig.txt'
        os.system(tempcomd)
        os.system("grep 'OW' waterout.pdbqt >> allwater.pdbqt")
        os.system("grep 'RESULT' waterout.pdbqt >> result.log")
        os.remove('waterout.pdbqt')

    os.remove('vinaconfig.txt')
    coods = np.genfromtxt('allwater.pdbqt', usecols=(5, 6, 7), dtype=float)
    energy = np.genfromtxt('result.log', usecols=3, dtype=float)

    selectcoods = np.compress(energy < -0.5, coods, axis=0)

    fit1 = scipy.cluster.hierarchy.fclusterdata(
        selectcoods, 0.5, criterion='distance', metric='euclidean')
    fit1 = fit1.astype(int)
    numclust1 = np.max(fit1)

    clustcoods = np.zeros((numclust1, 3), dtype=float)
    for i in range(1, numclust1+1):

        temp = np.compress(fit1 == i, selectcoods, axis=0)
        tempavg = np.mean(temp, axis=0)
        clustcoods[i-1, :] = tempavg

    fit2 = scipy.cluster.hierarchy.fclusterdata(
        clustcoods, 1.6, criterion='distance', metric='euclidean')
    fit2 = fit2.astype(int)
    numclust2 = np.max(fit2)

    finalcoods = np.zeros((numclust2, 3), dtype=float)
    for i in range(1, numclust2+1):

        temp = np.compress(fit2 == i, clustcoods, axis=0)
        tempavg = np.mean(temp, axis=0)
        finalcoods[i-1, :] = tempavg

    if os.path.isfile('predictedwaters.pdb'):
        os.rename('predictedwaters.pdb', '#predictedwaters.pdb')
    write_waterpdb('predictedwaters.pdb', finalcoods)

    os.remove('allwater.pdbqt')
    os.remove('result.log')
    os.remove('water.pdbqt')

    cmd.set("retain_order", 1)
    cmd.set("pdb_use_ter_records", 0)
    cmd.load(proteinfile, 'pro')
    cmd.show_as('cartoon', 'pro')
    cmd.color('blue', 'pro')
    cmd.load('predictedwaters.pdb', 'wats')
    cmd.color('red', 'wats')
    cmd.center('wats')


def option2():

    global dialog
    vinacomd = str(checkutilities())
    dialog = inputdata2()

    if dialog.result() == 1:
        runholowaterdock(vinacomd, dialog.proteinfile, dialog.ligandfile)


def runholowaterdock(vinacomd, proteinfile, ligandfile):

    from . import addwater
    addwater.main(ligandfile)

    waterfile()
    # Writes the water.pdbqt file

    from . import dockcheck
    dockcheck.main(proteinfile, ligandfile, vinacomd)

    os.remove('waterdetails.txt')
    os.remove('placedwaters.pdb')
    os.remove('water.pdbqt')

    cmd.set("retain_order", 1)
    cmd.set("pdb_use_ter_records", 0)
    cmd.load(proteinfile, 'pro')
    cmd.show_as('cartoon', 'pro')
    cmd.color('blue', 'pro')
    cmd.load('predictedwaters.pdb', 'wats')
    cmd.color('red', 'wats')
    cmd.load(ligandfile, 'lig')
    cmd.show_as('sticks', 'lig')
    cmd.center('lig')


def checkutilities():
    zz = distutils.spawn.find_executable('vina')
    homedir = str(os.path.expanduser('~'))
    file = os.path.join(homedir, 'pyvina.txt')

    if zz:  # Found vina under command 'vina'
        if os.path.isfile(zz):
            vinacomd = zz

    elif os.path.isfile(file):  # Previously stored the path to vina executible in pyvina.txt file
        fzz = open(file, 'r')
        vinacomd = fzz.read()
        fzz.close()

    else:  # Needs user to input path to vina executable

        vinapath_dialog = vinapath()
        vinapath_dialog.show()
        vinapath_dialog.raise_()

        if vinapath_dialog.result() == 1:
            fzz = open(file, 'r')
            vinacomd = fzz.read()
            fzz.close()

    return vinacomd


def waterfile():
    f1 = open('water.pdbqt', 'w')

    f1.write('REMARK  The pdbqt file for using water as a ligand\n')
    f1.write('ROOT\n')
    f1.write(
        'ATOM      1  OW  HOH   231       0.950  11.375  16.494  1.00  0.00    -0.411 OA\n')
    f1.write(
        'ATOM      2  HW1 HOH   231       1.766  11.375  17.071  1.00  0.00     0.205 HD\n')
    f1.write(
        'ATOM      3  HW2 HOH   231       0.134  11.375  17.071  1.00  0.00     0.205 HD\n')
    f1.write('ENDROOT\n')
    f1.write('TORSDOF 0')

    f1.close()


def write_waterpdb(waterfilename, coordinates):

    numatom = coordinates.shape[0]
    xyz = open(waterfilename, 'w')

    xyz.write("TITLE    Water Molecules Predicted by WaterDock\n")
    xyz.write("REMARK	Please Cite: Rapid and Accurate Prediction and Scoring of Water \
              Molecules in Protein Binding Sites\n")
    xyz.write("REMARK	DOI:10.1371/journal.pone.0032036\n")

    for i in range(0, numatom):
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
        xyz.write("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (
            header, serial, name, iCode, resname, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor))

    xyz.close()


class inputdata1(qtw.QDialog):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.setWindowTitle("Welcome to WaterDock 2.0")

        #  Interactive Widgets
        self.protein_file_input = qtw.QLineEdit()
        self.protein_file_choose_button = qtw.QPushButton('Choose')
        self.ligand_coordinates_inputted = qtw.QRadioButton('Enter center of box (A)')
        self.coordinates_widget = qtw.QWidget()
        self.xcom_input = qtw.QLineEdit()
        self.ycom_input = qtw.QLineEdit()
        self.zcom_input = qtw.QLineEdit()
        self.ligand_file_inputted = qtw.QRadioButton(
            'Coordinates from ligand file (mol2/pdb/pdbqt)')
        self.ligand_file_input = qtw.QLineEdit()
        self.ligand_file_choose_button = qtw.QPushButton('Choose')
        self.cancelrun_widget = qtw.QWidget()
        self.run_button = qtw.QPushButton('Run')
        self.cancel_button = qtw.QPushButton('Cancel')
        self.coordinates_widget.setEnabled(False)
        self.ligand_file_input.setEnabled(False)
        self.ligand_file_choose_button.setEnabled(False)
        self.ligopt = 0

        #  Layout
        layout = qtw.QGridLayout()
        layout.addWidget(qtw.QLabel('Predicting Waters in Apo-protein Structures'), 0, 0)
        layout.addWidget(qtw.QLabel('Import protein from File (pdbqt)'), 3, 0)
        layout.addWidget(self.protein_file_input, 3, 1)
        layout.addWidget(self.protein_file_choose_button, 3, 2)

        layout.addWidget(qtw.QLabel('Identify the ligand binding site'), 6, 0)
        self.coordinates_widget.setLayout(qtw.QHBoxLayout())
        self.coordinates_widget.layout().addWidget(qtw.QLabel('X'))
        self.coordinates_widget.layout().addWidget(self.xcom_input)
        self.coordinates_widget.layout().addWidget(qtw.QLabel('Y'))
        self.coordinates_widget.layout().addWidget(self.ycom_input)
        self.coordinates_widget.layout().addWidget(qtw.QLabel('Z'))
        self.coordinates_widget.layout().addWidget(self.zcom_input)
        layout.addWidget(self.ligand_coordinates_inputted, 7, 0)
        layout.addWidget(self.coordinates_widget, 7, 1)
        layout.addWidget(self.ligand_file_inputted, 8, 0)
        layout.addWidget(self.ligand_file_input, 8, 1)
        layout.addWidget(self.ligand_file_choose_button, 8, 2)

        self.cancelrun_widget.setLayout(qtw.QHBoxLayout())
        self.cancelrun_widget.layout().addWidget(self.run_button)
        self.cancelrun_widget.layout().addWidget(self.cancel_button)
        layout.addWidget(self.cancelrun_widget, 10, 0, 1, 4, qtc.Qt.AlignHCenter)

        layout.setSizeConstraint(layout.SetFixedSize)
        for row in range(0, 8):
            layout.setRowMinimumHeight(row, 15)
        self.setLayout(layout)

        # Connect signals to functions/slots
        self.protein_file_choose_button.clicked.connect(self.profilechoose)
        self.ligand_coordinates_inputted.clicked.connect(self.coordinpoptionchecked)
        self.ligand_file_inputted.clicked.connect(self.ligfileoptionchecked)
        self.ligand_file_choose_button.clicked.connect(self.ligfilechoose)
        self.run_button.clicked.connect(self.rungui)
        self.cancel_button.clicked.connect(self.reject)

        # Display
        self.exec()
        self.raise_()

    def profilechoose(self):
        self.protein_file_input.clear()
        file_filter = "pdbqt files (*.pdbqt)"
        profilenamechoose = qtw.QFileDialog().getOpenFileName(None, '', '', file_filter)[0]
        self.protein_file_input.insert(profilenamechoose)

    def ligfilechoose(self):
        self.ligand_file_input.clear()
        file_filter = "pdb files (*.pdb);; pdbqt files (*.pdbqt);; mol2 files (*.mol2)"
        ligfilenamechoose = qtw.QFileDialog().getOpenFileName(None, '', '', file_filter)[0]
        self.ligand_file_input.insert(ligfilenamechoose)

    def coordinpoptionchecked(self):
        self.coordinates_widget.setEnabled(True)
        self.ligand_file_input.setEnabled(False)
        self.ligand_file_choose_button.setEnabled(False)
        self.ligopt = int(1)

    def ligfileoptionchecked(self):
        self.coordinates_widget.setEnabled(False)
        self.ligand_file_input.setEnabled(True)
        self.ligand_file_choose_button.setEnabled(True)
        self.ligopt = int(2)

    def rungui(self):

        self.proteinfile = ''
        self.centerx = None
        self.centery = None
        self.centerz = None

        if os.path.isfile(self.protein_file_input.text()):
            self.proteinfile = self.protein_file_input.text()
        else:
            warning_proteinfile = qtw.QMessageBox(2, 'Missing File', 'Absolute path to protein \
                                                  File has not been inputted or file does not exist')
            warning_proteinfile.exec()

        if self.ligopt == 1:
            if self.xcom_input.text() and self.ycom_input.text() and self.zcom_input.text():
                try:
                    self.centerx = float(self.xcom_input.text())
                    self.centery = float(self.ycom_input.text())
                    self.centerz = float(self.zcom_input.text())
                except ValueError:
                    warning_coord_NaN = qtw.QMessageBox(2, 'Incorrect Coordinate(s)', 'Please use \
                                                        integer or decimal for every ligand binding site coordinate')
                    warning_coord_NaN.exec()
            else:
                warning_missing_coord = qtw.QMessageBox(2, 'Missing Coordinate(s)',
                                                        'Please specify X, Y, Z for center of box')
                warning_missing_coord.exec()

        elif self.ligopt == 2:
            if os.path.isfile(self.ligand_file_input.text()):
                ligandfile = self.ligand_file_input.text()

                UL = MDAnalysis.Universe(ligandfile)

                if ligandfile[-5:] == 'pdbqt':
                    heavylig = UL.select_atoms('not type HD')
                else:
                    heavylig = UL.select_atoms('not type H')

                com = np.mean(heavylig.positions, axis=0)
                self.centerx = float(com[0])
                self.centery = float(com[1])
                self.centerz = float(com[2])

            else:
                warning_ligandfile = qtw.QMessageBox(2, 'Missing File', 'Absolute path to \
                                                     ligand file has not been inputted or file does not exist')
                warning_ligandfile.exec()

        else:
            warning_ligoption = qtw.QMessageBox(2, 'Missing Ligand Binding Site Information',
                                                'Please select an option for identifying the ligand binding site')
            warning_ligoption.exec()

        coord_test = isinstance(self.centerx, float) and isinstance(self.centery, float) \
            and isinstance(self.centerz, float)
        if os.path.isfile(self.proteinfile) and coord_test:
            self.accept()


class inputdata2(qtw.QDialog):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.setWindowTitle("Welcome to WaterDock 2.0")

        # Interactive widgets
        self.protein_file_input = qtw.QLineEdit()
        self.protein_file_choose_button = qtw.QPushButton('Choose')
        self.ligand_file_input = qtw.QLineEdit()
        self.ligand_file_choose_button = qtw.QPushButton('Choose')
        self.run_button = qtw.QPushButton('Run')
        self.cancel_button = qtw.QPushButton('Cancel')

        # Layout
        layout = qtw.QGridLayout()
        layout.addWidget(qtw.QLabel("Predicting Waters in Holo-protein Structures"), 0, 0, 1, 5)
        layout.addWidget(qtw.QLabel("Import protein from File (pdbqt)"), 1, 0)
        layout.addWidget(self.protein_file_input, 1, 1)
        layout.addWidget(self.protein_file_choose_button, 1, 2)
        layout.addWidget(qtw.QLabel("Import ligand from File (pdb/pdbqt/mol2)"), 2, 0)
        layout.addWidget(self.ligand_file_input, 2, 1)
        layout.addWidget(self.ligand_file_choose_button, 2, 2)
        layout.addWidget(self.run_button, 3, 3)
        layout.addWidget(self.cancel_button, 3, 4)
        layout.setSizeConstraint(layout.SetFixedSize)
        self.setLayout(layout)

        # Connect signals to functions/slots
        self.protein_file_choose_button.clicked.connect(self.profilechoose)
        self.ligand_file_choose_button.clicked.connect(self.ligfilechoose)
        self.run_button.clicked.connect(self.rungui)
        self.cancel_button.clicked.connect(self.reject)

        # Display
        self.exec()
        self.raise_()

    def profilechoose(self):
        self.protein_file_input.clear()
        file_filter = "pdbqt files (*.pdbqt)"  # required input for Vina
        protfilenamechoose = qtw.QFileDialog().getOpenFileName(None, '', '', file_filter)[0]
        self.protein_file_input.insert(protfilenamechoose)

    def ligfilechoose(self):
        self.ligand_file_input.clear()
        file_filter = "pdb files (*.pdb);; pdbqt files (*.pdbqt);; mol2 files (*.mol2)"
        ligfilenamechoose = qtw.QFileDialog.getOpenFileName(None, '', '', file_filter)[0]
        self.ligand_file_input.insert(ligfilenamechoose)

    def rungui(self):

        self.proteinfile = ''
        self.ligandfile = ''

        if os.path.isfile(self.protein_file_input.text()):
            self.proteinfile = self.protein_file_input.text()
        else:
            warning_proteinfile = qtw.QMessageBox(2, 'Missing File',
                                                  'Absolute path to protein file has \
                                                  not been inputted or file does not exist')
            warning_proteinfile.exec()

        if os.path.isfile(self.ligand_file_input.text()):
            self.ligandfile = self.ligand_file_input.text()
        else:
            warning_ligandfile = qtw.QMessageBox(2, 'Missing File',
                                                 'Absolute path to ligand file has not been \
                                                 inputted or file does not exist')
            warning_ligandfile.exec()

        if os.path.isfile(self.proteinfile) and os.path.isfile(self.ligandfile):
            self.accept()


class vinapath(qtw.QDialog):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.setWindowTitle('Unable to find Vina executable')

        # Widgets
        self.vinapath_input = qtw.QLineEdit()
        self.choose_button = qtw.QPushButton('Choose')
        self.cancel_button = qtw.QPushButton('Cancel')
        self.okay_button = qtw.QPushButton('OK')

        # Layout
        layout = qtw.QGridLayout()
        layout.addWidget(qtw.QLabel("Enter absolute path to the vina executable."
                                    + "\nPath will be written to file 'pyvina.txt'"
                                    + "to save for future runs"), 0, 0, 1, 4)
        layout.addWidget(self.vinapath_input, 1, 0)
        layout.addWidget(self.choose_button, 1, 1)
        layout.addWidget(self.cancel_button, 2, 2)
        layout.addWidget(self.okay_button, 2, 3)
        layout.setSizeConstraint(layout.SetFixedSize)
        self.setLayout(layout)

        #  Connect signals to functions
        self.choose_button.clicked.connect(self.vinafilechoose)
        self.cancel_button.clicked.connect(self.reject)
        self.okay_button.clicked.connect(self.okay)

        self.exec()
        self.raise_()

    def vinafilechoose(self):
        self.vinapath_input.clear()
        vinafilenamechoose = qtw.QFileDialog().getOpenFileName()[0]
        self.vinapath_input.insert(vinafilenamechoose)

    def okay(self):

        if os.path.isfile(self.vinapath_input.text()):
            self.vinapathway = self.vinapath_input.text()
            homedir = str(os.path.expanduser('~'))
            file = os.path.join(homedir, 'pyvina.txt')
            f1 = open(file, 'w')
            f1.write(str(self.vinapathway))
            f1.close()
            self.accept()
        else:
            warning = qtw.QMessageBox(2, 'Missing Path', 'Path to vina executable has not been specified')
            warning.exec()
