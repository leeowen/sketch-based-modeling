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
        self.saveAsDat_button=QtWidgets.QPushButton('Save as .dat')
        self.saveAsImg_button=QtWidgets.QPushButton('Save as image')
        self.close_button=QtWidgets.QPushButton('Close')
        
        self.standardEllipse_checkBox=QtWidgets.QCheckBox('standard ellipse')
        self.standardEllipse_checkBox.setChecked(False)
       
        self.generalizedEllipse_checkBox=QtWidgets.QCheckBox('generalized ellipse')
        self.generalizedEllipse_checkBox.setChecked(False)
        
        self.originalEllipse_checkBox=QtWidgets.QCheckBox('original ellipse')
        self.originalEllipse_checkBox.setChecked(False)
        
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
        self.file_button=QtWidgets.QPushButton()
        self.file_button.setIcon(QtGui.QIcon(':fileOpen.png'))# ':' tells Qt that the following is resource file path
        self.file_button.setToolTip("Select File")
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
        checkBox_layout.addWidget(self.originalEllipse_checkBox)
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
        button_layout.addWidget(self.saveAsDat_button)
        button_layout.addWidget(self.saveAsImg_button)
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
        self.originalEllipse_checkBox.toggled.connect(self.drawMode_originalEllipse)
        
    
    def update_visibility_originalEllipse_mode(self,checked):
        self.J_spinBox.setVisible(not checked)
        self.Ea_lineEdit.setVisible(not checked)
        self.Em_lineEdit.setVisible(not checked)
        self.segment_comboBox.setVisible(not checked)
        self.saveAsDat_button.setVisible(not checked)
    
    def show_file_selected_dialog(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        FILE_FILTERS="data(*.dat);;All Files(*.*)"
        selected_filter="data(*.dat)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_path,selected_filter=QtWidgets.QFileDialog.getOpenFileName(self, 'Select File',dirPath+'/data',FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            self.filePath_lineEdit.setText(file_path)
            self.readData(file_path)
            
            
    def readData(self,file_path):
        f=open(file_lineEdit,'r')
        f1=f.readlines()
        self.vertices=[]
        i=0
        for v in f1:
            pos=v.split()
            pos[0]=float(pos[0])
            pos[1]=float(pos[1])
            pos[2]=float(pos[2])
            self.vertices.append(pos)
            i+=1
        
        self.numVertices=i
        f.close()
        
        
    def closeFn(self):
        self.close()
        
          
    def createEllipse(self):
        if self.generalizedEllipse_checkBox.isChecked():
            self.canvas.setGeneralizedEllipseFlag(True)
        else:
            self.canvas.setGeneralizedEllipseFlag(False)                
            
        if self.originalEllipse_checkBox.isChecked():
            self.canvas.setOriginalizedEllipseFlag(True)
        else:
            self.canvas.setOriginalizedEllipseFlag(False)
           
        if self.standardEllipse_checkBox.isChecked():
            self.canvas.setOriginalizedEllipseFlag(True)
        else:
            self.canvas.setOriginalizedEllipseFlag(False)
                
        
    def drawMode_standardEllipse(self,checked):
        """
        if cmds.control('standardEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('standardEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
        """
        if checked==True:
            pass
            #self.calculate_standardEllipse()
        else:
            # erase standard ellipse
            pass
                            
            
    def drawMode_generalizedEllipse(self,checked):
        """
        if cmds.control('generalizedEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('generalizedEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
        """

        self.J_spinBox.setVisible(checked)
        self.Ea_lineEdit.setVisible(checked)
        self.Em_lineEdit.setVisible(checked)
        self.segment_comboBox.setVisible(checked)
        self.saveAsDat_button.setVisible(checked)           
           
        
    def drawMode_originalEllipse(self,checked):
        self.canvas.draw_originalEllipse=checked
        # verify if generalised ellipse is also required to be draw
        # if yes, then no need to call update_visibiblity_originalEllipse_mode()
        flag=not self.generalizedEllipse_checkBox.isChecked()
        if checked and flag:
            self.update_visibility_originalEllipse_mode(checked)
            
              


class Canvas(QtWidgets.QDialog):
    backgroundColor=QtCore.Qt.white
    def __init__(self,parent=None):
        super(Canvas,self).__init__(parent)
        self.setMinimumSize(600,600)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
    def paintEvent(self,evt):
        painter=QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        painter.setBrush(Canvas.backgroundColor)
        painter.drawRect(self.rect())
        #self.calculate_generalizedEllipse()
        #self.calculate_originalEllipse()
        #self.calculate_standardEllipse()
        
        
    def calculate_standardEllipse(self):
        penColor=QtCore.Qt.Blue   
                  

    def calculate_generalizedEllipse(self):
        penColor=QtCore.Qt.Green
        
    
    def calculate_originalEllipse(self):
        penColor=QtCore.Qt.Red
           

if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("generalizedEllipseWin",exists=True):
        cmds.deleteUI("generalizedEllipseWin",wnd=True)
        
    w=CurveFittingWindowUI()
    w.show()  
    #w.raise_() 

