import maya.api.OpenMaya as om
import maya.cmds as cmds
from PySide2 import QtCore
from PySide2 import QtWidgets
from shiboken2 import wrapInstance
import maya.OpenMayaUI as omui
import math,sys

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy


def maya_useNewAPI():
    pass


def maya_main_window():
    main_window_ptr=omui.MQtUtil.mainWindow()
    return wrapInstance(long(main_window_ptr),QtWidgets.QWidget)
    
       
class WindowUI(QtWidgets.QMainWindow):
    
    def __init__(self,parent=maya_main_window()):
        super(WindowUI,self).__init__(parent)
        
        self.setWindowTitle("Generalized Ellipse")
        self.setMinimumWidth(800)
        self.setMinimumHeight(600)
        self.setObjectName("generalizedEllipseWin")
        
        self.create_widgets()
        self.create_layout()
        self.create_connections() 
        
        
    def create_widgets(self):
        self.button_generate=QtWidgets.QPushButton('Generate')
        
    def create_layout(self):
        layout_widget=QtWidgets.QWidget()
        self.setCentralWidget(layout_widget)
        main_layout=QtWidgets.QHBoxLayout(layout_widget)
        main_layout.addWidget(self.button_generate)
        
    def create_connections(self):
        self.button_generate.clicked.connect(self.createEllipse)
    
    def createEllipse(self):
        print "ellipse created"

if __name__=="__main__":
    if cmds.window("generalizedEllipseWin",exists=True):
        cmds.deleteUI("generalizedEllipseWin",wnd=True)
        
    w=WindowUI()
    w.show()    
    
    
#cmds.flushUndo()
#cmds.scriptEditorInfo(clearHistory=True)