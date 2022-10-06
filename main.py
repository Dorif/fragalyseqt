from PyQt5 import QtWidgets
import fragalyseqt, sys
class FragalyseApp(QtWidgets.QMainWindow, fragalyseqt.Ui_MainWindow):
    def __init__(self, parent=None):
        super(FragalyseApp, self).__init__(parent)
        self.setupUi(self)
def main():
    app = QtWidgets.QApplication(sys.argv)
    form = FragalyseApp()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
