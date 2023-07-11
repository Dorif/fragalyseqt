try:
    from PyQt5.QtWidgets import QMessageBox
except ImportError:
    from PyQt6.QtWidgets import QMessageBox
def msgbox(headerstr, msgstr, msgtype):
    mbox = QMessageBox()
    if msgtype == 1:
        mbox.setIcon(QMessageBox.Warning)
    elif msgtype == 2:
        mbox.setIcon(QMessageBox.Critical)
    else:
        mbox.setIcon(QMessageBox.Information)
    mbox.setWindowTitle(headerstr)
    mbox.setText(msgstr)
#Yes, this button is shown by default and is deafult button by default, but better I'll do it in explicit manner.
    mbox.setStandardButtons(QMessageBox.Ok)
    mbox.setDefaultButton(QMessageBox.Ok)
    mbox.exec()
