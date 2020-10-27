#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-
import sys, os
import numpy as np 

from PyQt5 import QtWidgets, QtGui, QtCore
from GUI_intrinsic_rate_constants_PyQt5 import Ui_MainWindow
import RateConstantCalculator as RCC 

import sys

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('PLQY calculator')

        # Load Button
        self.ui.comboBox.activated.connect(self.comboBox_activated)

        # Print File Directory Button
        self.ui.textEdit_PF.textChanged.connect(self.textEdit_PF_textChanged)
        self.ui.textEdit_DF.textChanged.connect(self.textEdit_DF_textChanged)

        # Calculate Button
        self.ui.textEdit_Phi_PF.textChanged.connect(self.textEdit_Phi_PF_textChanged)
        self.ui.textEdit_Phi_DF.textChanged.connect(self.textEdit_Phi_DF_textChanged)

        # Save Button
        self.ui.pushButton_calculate.clicked.connect(self.pushButton_calculate_clicked)

        # data initialize
        self.tauPF = 10e-9 
        self.tauDF = 10e-6

        self.kPF = 1/self.tauPF
        self.kDF = 1/self.tauDF

        self.Phi_PF = 0.5
        self.Phi_DF = 0.5
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
    def comboBox_activated(self):
        Index = self.ui.comboBox.currentIndex()
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
    def textEdit_PF_textChanged(self):
        string = self.ui.textEdit_PF.toPlainText()
        if not self.IsFloat(string):
            print('Prompt fluorescence value should be a positive value.')
            return 
        value = float(string)
        Index = self.ui.comboBox.currentIndex()
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
        Index = self.ui.comboBox.currentIndex()
        if Index == 0:  # lifetime
            self.tauDF = value*1e-6
        else: # rate constant
            self.kDF = value*1e5
    def textEdit_Phi_PF_textChanged(self):
        string = self.ui.textEdit_Phi_PF.toPlainText()
        if not self.IsFloat(string):
            print('Prompt fluorescence value should be a positive value.')
            return 
        self.Phi_PF = float(string) * 0.01
    def textEdit_Phi_DF_textChanged(self):
        string = self.ui.textEdit_Phi_DF.toPlainText()
        if not self.IsFloat(string):
            print('Delayed fluorescence value should be a positive value.')
            return 
        self.Phi_DF = float(string) * 0.01    
    def pushButton_calculate_clicked(self):
        # print Information
        Index = self.ui.comboBox.currentIndex()
        if Index == 0: # lifetime
            print('Mode   : Lifetime')
            print('tauPF  : {0} ns'.format( round(self.tauPF*1e9, 5) ) )
            print('tauDF  : {0} us'.format( round(self.tauDF*1e6, 5) ) )
        else: # rate constant
            print('Mode   : Rate Constant')
            print('kPF    : {0} x10^8 (1/s)'.format( round(self.kPF/1e8, 5) ) )
            print('kDF    : {0} x10^5 (1/s)'.format( round(self.kDF/1e5, 5) ) )
        print('Phi PF : {0} %'.format(self.Phi_PF*100) )
        print('Phi DF : {0} %'.format(self.Phi_DF*100) )
        print('')
        
        # save file path
        initialdir, extension = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", os.getcwd(), "All Files (*);;Text Files (*.txt)" )
        fpath = os.path.dirname(initialdir)
        fname = initialdir[(len(fpath)+1)::]
        if len(fname)==0 or len(fpath)==0:
            return

        # run
        if Index == 0: # lifetime
            tau_PF, tau_DF = self.tauPF, self.tauDF
        else:
            tau_PF, tau_DF = RCC.k2tau([self.kPF, self.kDF])
        RCC.script(tau_PF, tau_DF, self.Phi_PF, self.Phi_DF, fpath=fpath, fname=fname)

if __name__ == '__main__':
     app = QtWidgets.QApplication([])
     window = MainWindow()
     window.show()
     sys.exit(app.exec_())












