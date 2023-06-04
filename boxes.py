import sys
from PyQt5 import QtWidgets
def msgbox(headerstr, msgstr, msgtype):
    mbox = QtWidgets.QMessageBox()
    if msgtype == 1:
        mbox.setIcon(QtWidgets.QMessageBox.Warning)
    elif msgtype == 2:
        mbox.setIcon(QtWidgets.QMessageBox.Critical)
    else:
        mbox.setIcon(QtWidgets.QMessageBox.Information)
    mbox.setWindowTitle(headerstr)
    mbox.setText(msgstr)
#Yes, this button is shown by default and is deafult button by default, but better I'll do it in explicit manner.
    mbox.setStandardButtons(QtWidgets.QMessageBox.Ok)
    mbox.setDefaultButton(QtWidgets.QMessageBox.Ok)
    mbox.exec_()
