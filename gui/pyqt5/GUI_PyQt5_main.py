#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-
import sys, os
import numpy as np 

from PyQt5 import QtWidgets, QtGui, QtCore
from GUI_intrinsic_rate_constants_PyQt5 import Ui_MainWindow


### my module
guipath = os.path.dirname(os.path.abspath(__file__))
basepath = os.path.dirname(guipath)
if basepath not in sys.path:
    sys.path.append(basepath)
from src import rate_constants_calculator as rcc 


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('Rate constant calculator')

        # Load Button
        self.ui.comboBox_LF_k.activated.connect(self.comboBox_LF_k_activated)
        self.ui.comboBox_QY_B.activated.connect(self.comboBox_QY_B_activated)

        # lifetime/rate constants
        self.ui.textEdit_PF.textChanged.connect(self.textEdit_PF_textChanged)
        self.ui.textEdit_DF.textChanged.connect(self.textEdit_DF_textChanged)

        # Q.Y./B
        self.ui.textEdit_Phi_PF.textChanged.connect(self.textEdit_Phi_PF_textChanged)
        self.ui.textEdit_Phi_DF.textChanged.connect(self.textEdit_Phi_DF_textChanged)

        self.ui.textEdit_PLQY.textChanged.connect(self.textEdit_PLQY_textChanged)
        self.ui.label_PLQY.hide()
        self.ui.textEdit_PLQY.hide()
        self.ui.label_PLQY_unit.hide()

        # Save Button
        self.ui.pushButton_calculate.clicked.connect(self.pushButton_calculate_clicked)

        # data initialize
        self.tauPF = 15e-9 
        self.tauDF = 2.9e-6

        self.kPF = 1/self.tauPF
        self.kDF = 1/self.tauDF

        self.Phi_PF = 0.70
        self.Phi_DF = 0.12

        self.B_PF = 1.0
        self.B_DF = 1.0
        self.PLQY = 1.0
    
    # help functions
    def IsInteger(self,strr):
        try:
            datasize = int(strr)
            if datasize<1:
                return False
            return True
        except:
            return False
    
    def IsFloat(self,strr):
        try:
            datasize = float(strr)
            return True
        except:
            return False
    
    # clicked function
    def comboBox_LF_k_activated(self):
        Index = self.ui.comboBox_LF_k.currentIndex()
        if Index == 0: # lifetime
            self.ui.label_PF_unit.setText('ns')
            self.ui.label_DF_unit.setText('us')

            self.ui.textEdit_PF.setPlainText( str( round(self.tauPF*1e9, 5) ) )
            self.ui.textEdit_DF.setPlainText( str( round(self.tauDF*1e6, 5) ) )
        else: # rate constant
            self.ui.label_PF_unit.setText('x10^8 (1/s)')
            self.ui.label_DF_unit.setText('x10^5 (1/s)')

            self.ui.textEdit_PF.setPlainText( str( round(self.kPF/1e8, 5) ) )
            self.ui.textEdit_DF.setPlainText( str( round(self.kDF/1e5, 5) ) )
    
    def comboBox_QY_B_activated(self):
        Index = self.ui.comboBox_QY_B.currentIndex()
        if Index == 0: # Q.Y.
            self.ui.label_PF_unit_2.setText('%')
            self.ui.label_DF_unit_2.setText('%')

            self.ui.textEdit_Phi_PF.setPlainText( str( round(self.Phi_PF*100, 5) ) )
            self.ui.textEdit_Phi_DF.setPlainText( str( round(self.Phi_DF*100, 5) ) )
            
            self.ui.label_PLQY.hide()
            self.ui.textEdit_PLQY.hide()
            self.ui.label_PLQY_unit.hide()
        else: # B
            self.ui.label_PF_unit_2.setText('')
            self.ui.label_DF_unit_2.setText('')

            self.ui.textEdit_Phi_PF.setPlainText( str( round(self.B_PF, 5) ) )
            self.ui.textEdit_Phi_DF.setPlainText( str( round(self.B_DF, 5) ) )
            
            self.ui.label_PLQY.show()
            self.ui.textEdit_PLQY.show()
            self.ui.label_PLQY_unit.show()
   
    def textEdit_PF_textChanged(self):
        string = self.ui.textEdit_PF.toPlainText()
        if not self.IsFloat(string):
            print('Prompt fluorescence value should be a positive value.')
            return 
        value = float(string)
        Index = self.ui.comboBox_LF_k.currentIndex()
        if Index == 0:  # lifetime
            self.tauPF = value*1e-9
        else: # rate constant
            self.kPF = value*1e8
    
    def textEdit_DF_textChanged(self):
        string = self.ui.textEdit_DF.toPlainText()
        if not self.IsFloat(string):
            print('Delayed fluorescence value should be a positive value.')
            return 
        value = float(string)
        Index = self.ui.comboBox_LF_k.currentIndex()
        if Index == 0:  # lifetime
            self.tauDF = value*1e-6
        else: # rate constant
            self.kDF = value*1e5
    
    def textEdit_Phi_PF_textChanged(self):
        string = self.ui.textEdit_Phi_PF.toPlainText()
        if not self.IsFloat(string):
            print('Prompt fluorescence value should be a positive value.')
            return 
        Index = self.ui.comboBox_QY_B.currentIndex()
        if Index == 0: # Q.Y.
            self.Phi_PF = float(string) * 0.01
        else: # B
           self.B_PF = float(string)
    
    def textEdit_Phi_DF_textChanged(self):
        string = self.ui.textEdit_Phi_DF.toPlainText()
        if not self.IsFloat(string):
            print('Delayed fluorescence value should be a positive value.')
            return 
        Index = self.ui.comboBox_QY_B.currentIndex()
        if Index == 0: # Q.Y.
            self.Phi_DF = float(string) * 0.01  
        else: # B
           self.B_DF = float(string)  
    
    def textEdit_PLQY_textChanged(self):
        string = self.ui.textEdit_PLQY.toPlainText()
        if not self.IsFloat(string):
            print('PLQY should be a positive value.')
            return 
        self.PLQY = float(string) * 0.01
    
    def pushButton_calculate_clicked(self):
        # print Information
        Index1 = self.ui.comboBox_LF_k.currentIndex()
        if Index1 == 0: # lifetime
            print('Mode   : Lifetime')
            print('tauPF  : {0} ns'.format( round(self.tauPF*1e9, 5) ) )
            print('tauDF  : {0} us'.format( round(self.tauDF*1e6, 5) ) )
        else: # rate constant
            print('Mode   : Rate Constant')
            print('kPF    : {0} x10^8 (1/s)'.format( round(self.kPF/1e8, 5) ) )
            print('kDF    : {0} x10^5 (1/s)'.format( round(self.kDF/1e5, 5) ) )
        Index2 = self.ui.comboBox_QY_B.currentIndex()
        if Index2 == 0: # Q.Y.
            print('Mode   : Quantum Yield')
            print('Phi PF : {0} %'.format(self.Phi_PF*100) )
            print('Phi DF : {0} %'.format(self.Phi_DF*100) )
            print('')
        else:
            print('Mode   : B & PLQY')
            print('PLQY : {0} %'.format(self.PLQY*100) )
            print('B PF : {0}'.format(self.B_PF) )
            print('B DF : {0}'.format(self.B_DF) )
            print('')
        
        # save file path
        initialdir, extension = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", os.getcwd(), "All Files (*)" )
        fpath = os.path.dirname(initialdir)
        fname = initialdir[(len(fpath)+1)::]
        if len(fname)==0 or len(fpath)==0:
            return

        # run
        if Index1 == 0: # lifetime
            tau_PF, tau_DF = self.tauPF, self.tauDF
        else:
            tau_PF, tau_DF = rcc.k2tau([self.kPF, self.kDF])
        if Index2 == 0: # Q.Y.
            Phi_PF, Phi_DF = self.Phi_PF, self.Phi_DF
        else:
            Phi_PF, Phi_DF = rcc.cal_phi_PF_DF(self.PLQY, tau_PF, tau_DF, self.B_PF, self.B_DF)

        rcc.script(tau_PF, tau_DF, Phi_PF, Phi_DF, fpath=fpath, fname=fname)

if __name__ == '__main__':
     app = QtWidgets.QApplication([])
     window = MainWindow()
     window.show()
     sys.exit(app.exec_())












