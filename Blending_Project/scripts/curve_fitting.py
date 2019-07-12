import maya.api.OpenMaya as om
import maya.cmds as cmds
from PySide2 import QtCore, QtGui, QtWidgets
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
    
             
class CurveFittingWindowUI(QtWidgets.QMainWindow):
    
    def __init__(self,parent=maya_main_window()):
        super(CurveFittingWindowUI,self).__init__(parent)
        
        self.setWindowTitle("Generalized Ellipse")
        self.setMinimumWidth(800)
        self.setMinimumHeight(600)
        self.resize(800,600)
        self.setObjectName("generalizedEllipseWin")
        self.setWindowFlags(QtCore.Qt.Tool)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
               
        self.create_widgets()   
        self.create_layout()
        self.create_connection()
        
        #self.canvas.installEventFilter()
            
        
    def create_widgets(self):
        self.button_generate=QtWidgets.QPushButton('Generate')
        self.button_generate.setObjectName('generate_button')
        self.button_save=QtWidgets.QPushButton('Save')
        self.button_save.setObjectName('save_button')
        
        self.standardEllipse_checkBox=QtWidgets.QCheckBox('standard ellipse')
        self.standardEllipse_checkBox.setChecked(False)
        self.standardEllipse_checkBox.setMinimumWidth(200)
        self.standardEllipse_checkBox.setMaximumWidth(200)
        self.standardEllipse_checkBox.setObjectName('standardEllipse_checkBox')
        
        self.generalizedEllipse_checkBox=QtWidgets.QCheckBox('generalized ellipse')
        self.generalizedEllipse_checkBox.setChecked(False)
        self.generalizedEllipse_checkBox.setMinimumWidth(200)
        self.generalizedEllipse_checkBox.setMaximumWidth(200)
        self.generalizedEllipse_checkBox.setObjectName('generalizedEllipse_checkBox')
               
        #self.label_standardEllipse=QtWidgets.QLabel(':')


    def create_layout(self):            
        widget_layout_central=QtWidgets.QWidget()
        self.setCentralWidget(widget_layout_central)
        main_layout=QtWidgets.QHBoxLayout(widget_layout_central)
                
        # create layout for canvas
        left_layout=QtWidgets.QVBoxLayout()
        self.canvas=Canvas(self)
        left_layout.addWidget(self.canvas)
        # create layout for J, Ea, Em
               
        # create layout for checkbox and buttons
        right_layout=QtWidgets.QVBoxLayout()
        right_layout.addWidget(self.standardEllipse_checkBox)
        right_layout.addWidget(self.generalizedEllipse_checkBox)
        
        # add a spacer
        self.spacer=QtWidgets.QSpacerItem(0,100)
        right_layout.addSpacerItem(self.spacer)
        
        right_layout.addWidget(self.button_generate)
        right_layout.addWidget(self.button_save)
        
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        

    def create_connection(self):
        self.button_generate.clicked.connect(self.createEllipse)
        self.standardEllipse_checkBox.clicked.connect(self.drawMode_standardEllipse)
        self.generalizedEllipse_checkBox.clicked.connect(self.drawMode_generalizedEllipse)
        
          
    def createEllipse(self):
        print "ellipse created"
    
        
    def drawMode_standardEllipse(self):
        if cmds.control('standardEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('standardEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
            if value==True:
                self.drawStandardEllipse()
            else:
                pass
            
            
    def drawMode_generalizedEllipse(self):
        if cmds.control('generalizedEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('generalizedEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
            if value==True:
                self.drawGeneralizedEllipse()
            else:
                pass
        
                
    def drawStandardEllipse(self):
        print "draw standard ellipse"             


    def drawGeneralizedEllipse(self):
        print "draw generalized ellipse"


class Canvas(QtWidgets.QWidget):
    backgroundColor=QtGui.QColor(100,122,200)
    def __init__(self,parent):
        super(Canvas,self).__init__(parent)
        self.show()
        
    def paintevent(self,evt):
        painter=QtGui.QPainter()
        painter.begin(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        painter.setBrush(Canvas.backgroundColor)
        painter.drawRect(self.rect())
        painter.end(self)
    
    """
    def eventFilter(self,qObj,evt):
        if evt.type==QtCore.QEvent.Paint:
            self.paintEvent(evt)
            return True
            
        return False
    """
            

if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("generalizedEllipseWin",exists=True):
        cmds.deleteUI("generalizedEllipseWin",wnd=True)
        
    w=CurveFittingWindowUI()
    w.show()  
    w.raise_() 


    
#cmds.flushUndo()
#cmds.scriptEditorInfo(clearHistory=True)