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
            kPF, kDF = RCC.tau2k([self.tauPF, self.tauDF])
        else:
            kPF, kDF = self.kPF, self.kDF

        # QY
        PLQY = self.Phi_PF+self.Phi_DF
        phi_Tnr_PL_Array = np.linspace(0, 1-PLQY, int(200) )

        # Array
        ks_Array, ksr_Array, ksnr_Array, kisc_Array, kt_Array, ktr_Array, ktnr_Array, krisc_Array = RCC.IntrinsicRateConstants(kPF, kDF, self.Phi_PF, self.Phi_DF, phi_Tnr_PL=phi_Tnr_PL_Array)

        phi_sr_Array, phi_snr_Array, phi_isc_Array = RCC.phi_sr_snr_isc(ksr_Array, ksnr_Array, kisc_Array)
        phi_tr_Array, phi_tnr_Array, phi_risc_Array= RCC.phi_tr_tnr_risc(ktr_Array, ktnr_Array, krisc_Array)

        print('=============================== {0} ==============================='.format(fname) )
        print('* Rate Constants [1/s] :')
        print('  kPF= {0:>5.2e}, kDF = {1:>5.2e}'.format(kPF, kDF))
        for idx in [0,-1]:
            print('** phi_Tnr_PL = {0:>6.2f}%'.format(phi_Tnr_PL_Array[idx]*100) )
            print('  ks = {0:>5.2e}, ksr = {1:>5.2e}, ksnr = {2:>5.2e}, kisc  = {3:>5.2e}'.format(ks_Array[idx], ksr_Array[idx], ksnr_Array[idx], kisc_Array[idx]) )
            print('  kt = {0:>5.2e}, ktr = {1:>5.2e}, ktnr = {2:>5.2e}, krisc = {3:>5.2e}'.format(kt_Array[idx], ktr_Array[idx], ktnr_Array[idx], krisc_Array[idx]) )
        print('')
        print('* Process Efficiency [%] :')
        for idx in [0,-1]:
            print('** phi_Tnr_PL = {0:>6.2f}%'.format(phi_Tnr_PL_Array[idx]*100) )
            print('  phi_sr = {0:>6.2f}, phi_snr = {1:>6.2f}, phi_isc = {2:>6.2f}'.format(phi_sr_Array[idx]*100, phi_snr_Array[idx]*100, phi_isc_Array[idx]*100) )
            print('  phi_tr = {0:>6.2f}, phi_tnr = {1:>6.2f}, phi_risc= {2:>6.2f}'.format(phi_tr_Array[idx]*100, phi_tnr_Array[idx]*100, phi_risc_Array[idx]*100) )
        print('')

        # write
        with open( os.path.join(fpath, fname+'.txt'), 'w') as file:
            file.write('{0:>15s} {1:>15s} {2:>15s} {3:>15s} {4:>15s} {5:>15s} {6:>15s} {7:>15s} {8:>15s} {9:>15s} {10:>15s} {11:>15s} {12:>15s} {13:>15s}\n'.format('ks(1/s)', 'ksr(1/s)', 'ksnr(1/s)', 'kisc(1/s)', 'kt(1/s)', 'ktr(1/s)', 'ktnr(1/s)', 'krisc(1/s)', 'phi_sr(%)', 'phi_snr(%)', 'phi_isc(%)', 'phi_tr(%)', 'phi_tnr(%)', 'phi_risc(%)') )
            for ii in range( len(phi_Tnr_PL_Array) ):
                file.write('{0:>15.5e} {1:>15.5e} {2:>15.5e} {3:>15.5e} {4:>15.5e} {5:>15.5e} {6:>15.5e} {7:>15.5e} {8:>15.5f} {9:>15.5f} {10:>15.5f} {11:>15.5f} {12:>15.5f} {13:>15.5f}\n'.format(ks_Array[ii], ksr_Array[ii], ksnr_Array[ii], kisc_Array[ii], kt_Array[ii], ktr_Array[ii], ktnr_Array[ii], krisc_Array[ii], phi_sr_Array[ii]*100, phi_snr_Array[ii]*100, phi_isc_Array[ii]*100, phi_tr_Array[ii]*100, phi_tnr_Array[ii]*100, phi_risc_Array[ii]*100) )
        
if __name__ == '__main__':
     app = QtWidgets.QApplication([])
     window = MainWindow()
     window.show()
     sys.exit(app.exec_())












