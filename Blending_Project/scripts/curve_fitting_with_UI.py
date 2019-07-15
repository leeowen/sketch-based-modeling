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
    
             
class CurveFittingWindowUI(QtWidgets.QWidget):
    
    def __init__(self,parent=maya_main_window()):
        super(CurveFittingWindowUI,self).__init__(parent)
        
        self.setWindowTitle("Generalized Ellipse")
        self.setMinimumSize(900,650)
        self.setObjectName("generalizedEllipseWin")
        self.setWindowFlags(QtCore.Qt.WindowType.Dialog)
                       
        self.create_widgets()   
        self.create_layout()
        self.create_connection()
                   
            
    def create_widgets(self):
        self.generate_button=QtWidgets.QPushButton('Generate')
        self.generate_button.setObjectName('generate_button')
        self.save_button=QtWidgets.QPushButton('Save')
        self.save_button.setObjectName('save_button')
        self.close_button=QtWidgets.QPushButton('Close')
        
        self.standardEllipse_checkBox=QtWidgets.QCheckBox('standard ellipse')
        self.standardEllipse_checkBox.setChecked(False)
        self.standardEllipse_checkBox.setFixedWidth(200)
        self.standardEllipse_checkBox.setObjectName('standardEllipse_checkBox')
        
        self.generalizedEllipse_checkBox=QtWidgets.QCheckBox('generalized ellipse')
        self.generalizedEllipse_checkBox.setChecked(False)
        self.generalizedEllipse_checkBox.setFixedWidth(200)
        self.generalizedEllipse_checkBox.setObjectName('generalizedEllipse_checkBox')
        
        self.J_spinBox=QtWidgets.QSpinBox()
        self.J_spinBox.setValue(10)
        self.J_spinBox.setFixedWidth(200)
        self.J_spinBox.setMinimum(1)
        self.J_spinBox.setSingleStep(1)
        
        self.Ea_lineEdit=QtWidgets.QLineEdit()
        self.Ea_lineEdit.setFixedWidth(200)

        self.Em_lineEdit=QtWidgets.QLineEdit()
        self.Em_lineEdit.setFixedWidth(200)
        
        self.canvas=Canvas()       
        
        self.file_label=QtWidgets.QLabel('File:')
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.file_label.setSizePolicy(sizePolicy)
        self.filePath_lineEdit=QtWidgets.QLineEdit()
        self.file_label.setBuddy(self.filePath_lineEdit)
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.filePath_lineEdit.setSizePolicy(sizePolicy)
        self.file_button=QtWidgets.QPushButton('...')
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        self.file_button.setSizePolicy(sizePolicy)
        
        self.segment_comboBox=QtWidgets.QComboBox()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.segment_comboBox.setSizePolicy(sizePolicy)
        self.segment_comboBox.addItem('single piece')
        self.segment_comboBox.addItem('segment')
        

    def create_layout(self):   
        """     
        widget_layout_central=QtWidgets.QWidget()
        self.setCentralWidget(widget_layout_central)
        main_layout=QtWidgets.QHBoxLayout(widget_layout_central)
        """
        main_layout=QtWidgets.QHBoxLayout(self)
                
        # create layout for canvas
        left_layout=QtWidgets.QVBoxLayout()
        
        left_layout.addWidget(self.canvas)
        left_layout.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)        
              
        # create layout for different parameters and operation buttons
        right_layout=QtWidgets.QVBoxLayout()
        #right_layout.addStretch()
        right_layout.setContentsMargins(2,2,3,3)
        #right_layout.setSpacing(10)
        
        file_layout=QtWidgets.QHBoxLayout()
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.filePath_lineEdit)
        file_layout.addWidget(self.file_button)
        right_layout.addLayout(file_layout)
        
        checkBox_layout=QtWidgets.QVBoxLayout()
        checkBox_layout.addWidget(self.standardEllipse_checkBox)
        checkBox_layout.addWidget(self.generalizedEllipse_checkBox)
        right_layout.addLayout(checkBox_layout)
        
        combo_layout=QtWidgets.QFormLayout()
        combo_layout.addRow('Mode:',self.segment_comboBox)
        right_layout.addLayout(combo_layout)
        
        form_layout=QtWidgets.QFormLayout()
        form_layout.addRow('J:',self.J_spinBox)
        form_layout.addRow('Ea:',self.Ea_lineEdit)
        form_layout.addRow('Em:',self.Em_lineEdit)    
                   
        right_layout.addLayout(form_layout)
        
        button_layout=QtWidgets.QVBoxLayout()
        button_layout.addWidget(self.generate_button)
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.close_button)
        right_layout.addLayout(button_layout)
        
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        

    def create_connection(self):
        self.file_button.clicked.connect(self.show_file_selected_dialog)
        self.generate_button.clicked.connect(self.createEllipse)
        self.close_button.clicked.connect(self.closeFn)
        self.standardEllipse_checkBox.toggled.connect(self.drawMode_standardEllipse)
        self.generalizedEllipse_checkBox.toggled.connect(self.drawMode_generalizedEllipse)
        
    def show_file_selected_dialog(self):
        pass
        
    def closeFn(self):
        self.close()
          
    def createEllipse(self):
        print "ellipse created"
    
        
    def drawMode_standardEllipse(self,checked):
        """
        if cmds.control('standardEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('standardEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
        """
        if checked==True:
            self.canvas.draw_standardEllipse=True
        else:
            # erase standard ellipse
            self.canvas.draw_standardEllipse=False
                            
            
    def drawMode_generalizedEllipse(self,checked):
        """
        if cmds.control('generalizedEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('generalizedEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
        """
        if checked==True:
            self.canvas.draw_generalizedEllipse=True
        else:
            # erase generalized ellipse
            self.canvas.draw_generalizedEllipse=False
        
                
    def drawStandardEllipse(self):
        print "draw standard ellipse"             


    def drawGeneralizedEllipse(self):
        print "draw generalized ellipse"


class Canvas(QtWidgets.QDialog):
    backgroundColor=QtCore.Qt.white
    def __init__(self,parent=None):
        super(Canvas,self).__init__(parent)
        self.setMinimumSize(600,600)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
        self.draw_standardEllipse=False
        self.draw_generalizedEllipse=False
        
    def paintEvent(self,evt):
        painter=QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        painter.setBrush(Canvas.backgroundColor)
        painter.drawRect(self.rect())
           

if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("generalizedEllipseWin",exists=True):
        cmds.deleteUI("generalizedEllipseWin",wnd=True)
        
    w=CurveFittingWindowUI()
    w.show()  
    #w.raise_() 

