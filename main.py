try:
    from PyQt5.QtWidgets import QMainWindow, QApplication
except ImportError:
    from PyQt6.QtWidgets import QMainWindow, QApplication
import fragalyseqt, sys
class FragalyseApp(QMainWindow, fragalyseqt.Ui_MainWindow):
    def __init__(self, parent=None):
        super(FragalyseApp, self).__init__(parent)
        self.setupUi(self)
def main():
    app = QApplication(sys.argv)
    form = FragalyseApp()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
