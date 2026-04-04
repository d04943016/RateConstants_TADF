#!/usr/bin/env python3
# Copyright (c) 2020 Wei-Kai Lee. All rights reserved

# coding=utf-8
# -*- coding: utf8 -*-
import sys
import os
import json
import numpy as np

# Resolve project root (4 levels up from src/ui/pyqt5/)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
# Add this directory for local imports (ui_window)
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from PyQt5 import QtWidgets, QtGui, QtCore
from ui_window import Ui_MainWindow
from src.engine import calculator as engine

# Load config
_CFG_PATH = os.path.join(THIS_DIR, 'config', 'config.json')
with open(_CFG_PATH) as _f:
    _CFG = json.load(_f)

DEFAULT_OUTPUT_DIR = os.path.join(ROOT, _CFG.get('default_output_dir', 'data/output'))
os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle(_CFG.get('window_title', 'Rate Constant Calculator'))

        # Connect signals
        self.ui.comboBox_LF_k.activated.connect(self.comboBox_LF_k_activated)
        self.ui.comboBox_QY_B.activated.connect(self.comboBox_QY_B_activated)

        self.ui.textEdit_PF.textChanged.connect(self.textEdit_PF_textChanged)
        self.ui.textEdit_DF.textChanged.connect(self.textEdit_DF_textChanged)

        self.ui.textEdit_Phi_PF.textChanged.connect(self.textEdit_Phi_PF_textChanged)
        self.ui.textEdit_Phi_DF.textChanged.connect(self.textEdit_Phi_DF_textChanged)

        self.ui.textEdit_PLQY.textChanged.connect(self.textEdit_PLQY_textChanged)
        self.ui.label_PLQY.hide()
        self.ui.textEdit_PLQY.hide()
        self.ui.label_PLQY_unit.hide()

        self.ui.pushButton_calculate.clicked.connect(self.pushButton_calculate_clicked)

        # Default values
        self.tauPF = 15e-9
        self.tauDF = 2.9e-6
        self.kPF = 1 / self.tauPF
        self.kDF = 1 / self.tauDF
        self.Phi_PF = 0.70
        self.Phi_DF = 0.12
        self.B_PF = 1.0
        self.B_DF = 1.0
        self.PLQY = 1.0

    def IsInteger(self, strr):
        try:
            datasize = int(strr)
            return datasize >= 1
        except:
            return False

    def IsFloat(self, strr):
        try:
            float(strr)
            return True
        except:
            return False

    def comboBox_LF_k_activated(self):
        Index = self.ui.comboBox_LF_k.currentIndex()
        if Index == 0:  # lifetime
            self.ui.label_PF_unit.setText('ns')
            self.ui.label_DF_unit.setText('us')
            self.ui.textEdit_PF.setPlainText(str(round(self.tauPF * 1e9, 5)))
            self.ui.textEdit_DF.setPlainText(str(round(self.tauDF * 1e6, 5)))
        else:           # rate constant
            self.ui.label_PF_unit.setText('x10^8 (1/s)')
            self.ui.label_DF_unit.setText('x10^5 (1/s)')
            self.ui.textEdit_PF.setPlainText(str(round(self.kPF / 1e8, 5)))
            self.ui.textEdit_DF.setPlainText(str(round(self.kDF / 1e5, 5)))

    def comboBox_QY_B_activated(self):
        Index = self.ui.comboBox_QY_B.currentIndex()
        if Index == 0:  # Q.Y.
            self.ui.label_PF_unit_2.setText('%')
            self.ui.label_DF_unit_2.setText('%')
            self.ui.textEdit_Phi_PF.setPlainText(str(round(self.Phi_PF * 100, 5)))
            self.ui.textEdit_Phi_DF.setPlainText(str(round(self.Phi_DF * 100, 5)))
            self.ui.label_PLQY.hide()
            self.ui.textEdit_PLQY.hide()
            self.ui.label_PLQY_unit.hide()
        else:           # B
            self.ui.label_PF_unit_2.setText('')
            self.ui.label_DF_unit_2.setText('')
            self.ui.textEdit_Phi_PF.setPlainText(str(round(self.B_PF, 5)))
            self.ui.textEdit_Phi_DF.setPlainText(str(round(self.B_DF, 5)))
            self.ui.label_PLQY.show()
            self.ui.textEdit_PLQY.show()
            self.ui.label_PLQY_unit.show()

    def textEdit_PF_textChanged(self):
        string = self.ui.textEdit_PF.toPlainText()
        if not self.IsFloat(string):
            return
        value = float(string)
        if self.ui.comboBox_LF_k.currentIndex() == 0:
            self.tauPF = value * 1e-9
        else:
            self.kPF = value * 1e8

    def textEdit_DF_textChanged(self):
        string = self.ui.textEdit_DF.toPlainText()
        if not self.IsFloat(string):
            return
        value = float(string)
        if self.ui.comboBox_LF_k.currentIndex() == 0:
            self.tauDF = value * 1e-6
        else:
            self.kDF = value * 1e5

    def textEdit_Phi_PF_textChanged(self):
        string = self.ui.textEdit_Phi_PF.toPlainText()
        if not self.IsFloat(string):
            return
        if self.ui.comboBox_QY_B.currentIndex() == 0:
            self.Phi_PF = float(string) * 0.01
        else:
            self.B_PF = float(string)

    def textEdit_Phi_DF_textChanged(self):
        string = self.ui.textEdit_Phi_DF.toPlainText()
        if not self.IsFloat(string):
            return
        if self.ui.comboBox_QY_B.currentIndex() == 0:
            self.Phi_DF = float(string) * 0.01
        else:
            self.B_DF = float(string)

    def textEdit_PLQY_textChanged(self):
        string = self.ui.textEdit_PLQY.toPlainText()
        if not self.IsFloat(string):
            return
        self.PLQY = float(string) * 0.01

    def pushButton_calculate_clicked(self):
        Index1 = self.ui.comboBox_LF_k.currentIndex()
        if Index1 == 0:
            print('Mode   : Lifetime')
            print('tauPF  : {} ns'.format(round(self.tauPF * 1e9, 5)))
            print('tauDF  : {} us'.format(round(self.tauDF * 1e6, 5)))
        else:
            print('Mode   : Rate Constant')
            print('kPF    : {} x10^8 (1/s)'.format(round(self.kPF / 1e8, 5)))
            print('kDF    : {} x10^5 (1/s)'.format(round(self.kDF / 1e5, 5)))

        Index2 = self.ui.comboBox_QY_B.currentIndex()
        if Index2 == 0:
            print('Mode   : Quantum Yield')
            print('Phi PF : {} %'.format(self.Phi_PF * 100))
            print('Phi DF : {} %'.format(self.Phi_DF * 100))
        else:
            print('Mode   : B & PLQY')
            print('PLQY : {} %'.format(self.PLQY * 100))
            print('B PF : {}'.format(self.B_PF))
            print('B DF : {}'.format(self.B_DF))
        print('')

        # Save file dialog - default to configured output dir
        initialdir, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save file", DEFAULT_OUTPUT_DIR, "All Files (*)"
        )
        if not initialdir:
            return
        fpath = os.path.dirname(initialdir)
        fname = os.path.basename(initialdir)
        if not fname or not fpath:
            return

        if Index1 == 0:
            tau_PF, tau_DF = self.tauPF, self.tauDF
        else:
            taus = engine.k2tau([self.kPF, self.kDF])
            tau_PF, tau_DF = float(taus[0]), float(taus[1])

        if Index2 == 0:
            Phi_PF, Phi_DF = self.Phi_PF, self.Phi_DF
        else:
            Phi_PF, Phi_DF = engine.cal_phi_PF_DF(self.PLQY, tau_PF, tau_DF, self.B_PF, self.B_DF)

        engine.script(tau_PF, tau_DF, Phi_PF, Phi_DF, fpath=fpath, fname=fname)


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
