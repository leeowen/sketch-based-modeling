import maya.api.OpenMaya as om
import maya.cmds as cmds
from PySide2 import QtCore, QtGui, QtWidgets
from shiboken2 import wrapInstance
import maya.OpenMayaUI as omui
import math,sys

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy as np


def maya_useNewAPI():
    pass


def maya_main_window():
    main_window_ptr=omui.MQtUtil.mainWindow()
    return wrapInstance(long(main_window_ptr),QtWidgets.QWidget)
    
             
class CurveFittingWindowUI(QtWidgets.QWidget):
    
    def __init__(self,parent=maya_main_window()):
        super(CurveFittingWindowUI,self).__init__(parent)
        
        self.setWindowTitle("Generalized Ellipse")
        self.setMinimumSize(900,700)
        self.setObjectName("generalizedEllipseWin")
        self.setWindowFlags(QtCore.Qt.WindowType.Dialog)
                       
        self.create_widgets()   
        self.create_layout()
        self.create_connection()
                                   
    
    def sizeHint(self):
        return QtCore.QSize(900,700)
        
        
    def minimumSize(self):
        return QtCore.Qsize(900,700)        
        
        
    def create_widgets(self):
        self.saveAsDat_button=QtWidgets.QPushButton('Save as .dat')
        self.saveAsImg_button=QtWidgets.QPushButton('Save as image')
        self.close_button=QtWidgets.QPushButton('Close')
        
        self.standardEllipse_checkBox=QtWidgets.QCheckBox('standard ellipse')
        self.standardEllipse_checkBox.setChecked(False)
        self.standardEllipse_checkBox.setVisible(False)
       
        self.generalizedEllipse_checkBox=QtWidgets.QCheckBox('generalized ellipse')
        self.generalizedEllipse_checkBox.setChecked(False)
        self.generalizedEllipse_checkBox.setVisible(False)
        
        self.originalEllipse_checkBox=QtWidgets.QCheckBox('original ellipse')
        self.originalEllipse_checkBox.setChecked(False)
        self.originalEllipse_checkBox.setVisible(False)
                       
        self.radio_group=QtWidgets.QGroupBox()  
        self.manualJ_mode_radioButton=QtWidgets.QRadioButton('manual J')
        self.autoJ_mode_radioButton=QtWidgets.QRadioButton('auto J') 
        
        self.J_label=QtWidgets.QLabel('J:')
        self.J_label.setFixedWidth(50)
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.J_label.setSizePolicy(sizePolicy)
        self.J_spinBox=QtWidgets.QSpinBox()
        self.J_spinBox.setValue(10)
        self.J_spinBox.setFixedWidth(100)
        self.J_spinBox.setMinimum(1)
        self.J_spinBox.setSingleStep(1)
        self.J_spinBox.setReadOnly(not self.manualJ_mode_radioButton.isChecked())

        self.Ea_label=QtWidgets.QLabel('Ea:') 
        self.Ea_label.setFixedWidth(50)   
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.Ea_label.setSizePolicy(sizePolicy)    
        self.Ea_lineEdit=QtWidgets.QLineEdit()
        self.Ea_lineEdit.setFixedWidth(100)
        self.Ea_lineEdit.setMaxLength(5)
        self.Ea_lineEdit.setReadOnly(True)

        self.Em_label=QtWidgets.QLabel('Em:')
        self.Em_label.setFixedWidth(50)
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.Em_label.setSizePolicy(sizePolicy)
        self.Em_lineEdit=QtWidgets.QLineEdit()
        self.Em_lineEdit.setFixedWidth(100)
        self.Em_lineEdit.setMaxLength(5)
        self.Ea_lineEdit.setReadOnly(True)
        
        self.canvas=Canvas()    
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        self.canvas.setSizePolicy(sizePolicy) 
        
        self.lineWidth_doubleSpinBox=QtWidgets.QDoubleSpinBox()
        self.lineWidth_doubleSpinBox.setRange(1.0,5.0)
        self.lineWidth_doubleSpinBox.setValue(1.0)
        self.lineWidth_doubleSpinBox.setSingleStep(0.5)  
        self.lineWidth_doubleSpinBox.setFixedWidth(120)
        
        self.file_label=QtWidgets.QLabel('File:')
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.file_label.setSizePolicy(sizePolicy)
        self.filePath_lineEdit=QtWidgets.QLineEdit()
        self.filePath_lineEdit.setMaximumWidth(150)
        self.file_label.setBuddy(self.filePath_lineEdit)
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.filePath_lineEdit.setSizePolicy(sizePolicy)
        self.file_button=QtWidgets.QPushButton()
        self.file_button.setIcon(QtGui.QIcon(':fileOpen.png'))# ':' tells Qt that the following is resource file path
        self.file_button.setToolTip("Select File")
        sizePolicy=QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        self.file_button.setSizePolicy(sizePolicy)
        
        self.segment_label=QtWidgets.QLabel('segment method:')
        self.segment_label.setMaximumWidth(150)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        self.segment_label.setSizePolicy(sizePolicy)
        self.segment_comboBox=QtWidgets.QComboBox()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        self.segment_comboBox.setSizePolicy(sizePolicy)
        self.segment_comboBox.setMaximumWidth(150)
        self.segment_comboBox.addItem('single piece')
        self.segment_comboBox.addItem('segment')
        
        self.radio_group.setVisible(False)
        self.segment_label.setVisible(False)
        self.segment_comboBox.setVisible(False)
        self.saveAsDat_button.setVisible(False)
        self.saveAsImg_button.setVisible(False)
        

    def create_layout(self):   
        main_layout=QtWidgets.QHBoxLayout(self)
                 
        center_layout=QtWidgets.QHBoxLayout()
        center_layout.addWidget(self.canvas)              
        center_layout.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)        
              
        # create layout for different parameters and operation buttons
        right_layout=QtWidgets.QVBoxLayout()
        
        right_layout.setContentsMargins(2,2,3,3)
        #right_layout.setSpacing(10)
        
        file_layout=QtWidgets.QHBoxLayout()
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.filePath_lineEdit)
        file_layout.addWidget(self.file_button)
        right_layout.addLayout(file_layout)
        
        lineWidth_layout=QtWidgets.QFormLayout()
        lineWidth_layout.addRow('line width:',self.lineWidth_doubleSpinBox)
        right_layout.addLayout(lineWidth_layout)
        
        right_layout.addStretch(1)
        
        checkBox_layout=QtWidgets.QVBoxLayout()
        checkBox_layout.addWidget(self.standardEllipse_checkBox)
        checkBox_layout.addWidget(self.generalizedEllipse_checkBox)
        checkBox_layout.addWidget(self.originalEllipse_checkBox)
        right_layout.addLayout(checkBox_layout)
        
        right_layout.addStretch(1)
        
        segmentMethod_layout=QtWidgets.QVBoxLayout()
        segmentMethod_layout.addWidget(self.segment_label)
        segmentMethod_layout.addWidget(self.segment_comboBox)
        right_layout.addLayout(segmentMethod_layout) 
        
        right_layout.addStretch(1)       
     
        grid_layout=QtWidgets.QGridLayout()
        grid_layout.addWidget(self.autoJ_mode_radioButton,0,0,1,2)
        grid_layout.addWidget(self.manualJ_mode_radioButton,1,0,1,2)
        grid_layout.addWidget(self.J_label,2,0)
        grid_layout.addWidget(self.J_spinBox,2,1)
        grid_layout.addWidget(self.Ea_label,3,0)
        grid_layout.addWidget(self.Ea_lineEdit,3,1)
        grid_layout.addWidget(self.Em_label,4,0)    
        grid_layout.addWidget(self.Em_lineEdit,4,1)
        self.radio_group.setLayout(grid_layout)
        
        right_layout.addWidget(self.radio_group)           
        
        right_layout.addStretch(1)
        
        button_layout=QtWidgets.QVBoxLayout()
        button_layout.addWidget(self.saveAsDat_button)
        button_layout.addWidget(self.saveAsImg_button)
        button_layout.addWidget(self.close_button)
        right_layout.addLayout(button_layout)
        
        center_layout.addLayout(right_layout)
        
        main_layout.addLayout(center_layout)
        

    def create_connection(self):
        self.saveAsDat_button.clicked.connect(self.save_dataFn)
        self.saveAsImg_button.clicked.connect(self.save_imageFn)
        self.file_button.clicked.connect(self.show_file_selected_dialog)
        self.close_button.clicked.connect(self.closeFn)
        self.standardEllipse_checkBox.toggled.connect(self.drawMode_standardEllipse)
        self.generalizedEllipse_checkBox.toggled.connect(self.drawMode_generalizedEllipse)
        self.originalEllipse_checkBox.toggled.connect(self.drawMode_originalEllipse)
        self.J_spinBox.valueChanged.connect(self.update_manual_J_value)
        self.lineWidth_doubleSpinBox.valueChanged.connect(self.canvas.setLineWidth)
        self.manualJ_mode_radioButton.toggled.connect(self.update_visibility_manual_J_mode)
        self.autoJ_mode_radioButton.toggled.connect(self.update_visibility_auto_J_mode)
        self.segment_comboBox.currentTextChanged.connect(self.segmentMode_change)
        
    
    def segmentMode_change(self,text):
        if text=='single piece':
            self.canvas.single_piece_mode=True
            self.canvas.segment_mode=False
        elif text=='segment':
            self.canvas.segment_mode=True
            self.canvas.single_piece_mode=False
            
        self.canvas.update()
                        
    
    def update_visibility_manual_J_mode(self,checked):
        self.canvas.manualJ_mode=checked
        self.canvas.autoJ_mode=not checked
        self.J_spinBox.setReadOnly(not checked)
        self.canvas.repaint()# update() prevents multiple fast repaints. call repaint() instead.
        self.showEaEm()
        
        
    def update_visibility_auto_J_mode(self,checked):
        self.canvas.autoJ_mode=checked
        self.canvas.manualJ_mode=not checked
        self.J_spinBox.setReadOnly(checked)
        self.canvas.repaint()
        self.showEaEm()
        self.J_spinBox.setValue(self.canvas.autoJ)
        
    
    def update_manual_J_value(self,J):
        self.canvas.setManualJ(J)
        self.showEaEm()
        
            
    def save_dataFn(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        dirPath+='data/'
        FILE_FILTERS="DAT(*.dat);All Files(*.*)"
        selected_filter="DAT(*.dat)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_path,selected_filter=QtWidgets.QFileDialog.getSaveFileName(self, 'save',dirPath+'data',FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            f=open(file_path,"w+")
            f.write('range:')
            f.write('0-360 \n')
            f.write('a: ')
            for i in range(self.canvas.numPt):
                f.write(str(self.a[i])+' ')
            f.write('\n')
            f.write('b: ')
            for i in range(self.canvas.numPt):
                f.write(str(self.b[i])+' ')
        f.close()
        
        
    def save_imageFn(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        dirPath+='images/'
        FILE_FILTERS="PNG(*.png);All Files(*.*)"
        selected_filter="PNG(*.png)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_path,selected_filter=QtWidgets.QFileDialog.getSaveFileName(self, 'save',dirPath+'images',FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            self.saveFile(file_path)
            
         
    def saveFile(self,file_path):
        file=QtCore.QFile(file_path)
        if file.open(QtCore.QIODevice.WriteOnly):
            pixmap=QtGui.QPixmap(self.canvas.size())
            self.canvas.render(pixmap)
            pixmap.save(file,"PNG")
            file.close()
        else:
            QMessageBox.warning(self,'curve fitting',tr("Cannot write file %1. \nError:%2").arg(file_path).arg(file.errorString()))
        
             
    def update_visibility_originalEllipse_mode(self,checked):
        self.radio_group.setVisible(not checked)
        self.segment_label.setVisible(not checked)
        self.segment_comboBox.setVisible(not checked)
        self.saveAsDat_button.setVisible(not checked)
        self.saveAsImg_button.setVisible(checked)


    def update_visibility_standardEllipse_mode(self,checked):
        self.saveAsImg_button.setVisible(not checked)
        
    
    def show_file_selected_dialog(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        FILE_FILTERS="data(*.dat);;All Files(*.*)"
        selected_filter="data(*.dat)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_path,selected_filter=QtWidgets.QFileDialog.getOpenFileName(self, 'Select File',dirPath+'data',FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            self.filePath_lineEdit.setText(file_path)
            self.standardEllipse_checkBox.setVisible(True)
            self.generalizedEllipse_checkBox.setVisible(True)
            self.originalEllipse_checkBox.setVisible(True)
            self.canvas.readFile(file_path)
            self.canvas.update()
            
        
    def closeFn(self):
        self.close() 
                
        
    def drawMode_standardEllipse(self,checked):
        self.canvas.standardEllipse=checked
        self.canvas.update()
        
        self.saveAsDat_button.setVisible(checked)
        self.saveAsImg_button.setVisible(checked)
                            
            
    def drawMode_generalizedEllipse(self,checked):
        """
        if cmds.control('generalizedEllipse'+'_checkBox',exists=True):
            ptr=omui.MQtUtil.findControl('generalizedEllipse_checkBox')
            checkBox=wrapInstance(long(ptr),QtWidgets.QCheckBox)
            value=checkBox.isChecked()
        """
        self.radio_group.setVisible(checked)
        self.segment_label.setVisible(checked)
        self.segment_comboBox.setVisible(checked)
        self.saveAsDat_button.setVisible(checked)  
        self.saveAsImg_button.setVisible(checked)  
        self.canvas.generalisedEllipse=checked
        self.canvas.update()
        
       
    def showEaEm(self):   
        s1=str(self.canvas.getEa()*100)
        s1=s1[0:4]
        s2=str(self.canvas.getEm()*100)
        s2=s2[0:4]
        self.Ea_lineEdit.setText(s1+'%')     
        self.Em_lineEdit.setText(s2+'%')
                  
        
    def drawMode_originalEllipse(self,checked):
        self.canvas.originalEllipse=checked
        self.canvas.update()

        self.saveAsDat_button.setVisible(checked)
        self.saveAsImg_button.setVisible(checked)
            
              
class Canvas(QtWidgets.QDialog):
    backgroundColor=QtCore.Qt.white
    Ea_criteria=0.01
    Em_criteria=0.05
    
    def __init__(self,parent=None):
        super(Canvas,self).__init__(parent)
        self.setMinimumSize(650,600)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.originalEllipse=False
        self.generalisedEllipse=False
        self.standardEllipse=False
        self.vertices=[]
        self.angles=[]
        self.d=[]
        self.d_bar=[]
        self.Ea=0.0
        self.Em=0.0
        self.numPt=0
        self.manualJ=10
        self.autoJ=0
        self.a=[]
        self.b=[]
        self.generalisedEllipseVertices=[]
        self.lineWidth=1
        self.manualJ_mode=False
        self.autoJ_mode=False
        self.single_piece_mode=True
        self.segment_mode=False
        
        
    def sizeHint(self):
        return QtCore.QSize(650,600)
        
        
    def minimumSize(self):
        return QtCore.Qsize(650,600)
        
        
    def setManualJ(self,j):
        self.manualJ=j
        self.update()
        
    
    def readFile(self,file_path):
        f=open(file_path,'r')
        content=f.readlines()
        self.center=QtCore.QPointF(0.0,0.0)
        self.numPt=0
        self.vertices=[]
        self.angles=[]
        self.d=[]
        self.d_bar=[]
        self.Ea=0.0
        self.Em=0.0
        for line in content:
            p=line.split()
            self.numPt+=1
            self.vertices.append(om.MVector(float(p[0])*300+self.width()/2.,float(p[1])*300,float(p[2])*300+self.height()/2.))
            self.center+=QtCore.QPointF(float(p[0])*300+self.width()/2.0,float(p[2])*300+self.height()/2.0)
        
        self.center=self.center/self.numPt
        for i in range(self.numPt):
            self.d_bar.append(math.sqrt((self.vertices[i][0]-self.center.x())**2+(self.vertices[i][2]-self.center.y())**2))
        self.getAngle()
        f.close()
        
        # split data for 2 segments respectively
        self.a_first_half=[]
        self.b_first_half=[]
        I=self.numPt/2
        self.vertices_first_half=self.vertices[0:I]
        self.angles_first_half=self.angles[0:I]
        self.vertices_second_half=self.vertices[I:self.numPt]
        self.vertices_second_half.append(self.vertices[0])
        self.angles_second_half=self.angles[I:self.numPt]
        self.angles_second_half.append(self.angles[0])
        self.center_first_half=QtCore.QPointF(0.0,0.0)
        self.center_second_half=QtCore.QPointF(0.0,0.0)
        for i in range(I):
            self.center_first_half+=QtCore.QPointF(self.vertices_first_half[i][0],self.vertices_first_half[i][2])
            self.center_second_half+=QtCore.QPointF(self.vertices_second_half[i][0],self.vertices_second_half[i][2])
        self.first_half_center/=(I+1)
        self.second_half_center/=(I+1)
            
    
    def getAngle(self):
        for i in range(0,self.numPt):
            anglem=(self.vertices[i][2]-self.center.y())/math.sqrt((self.vertices[i][2]-self.center.y())**2+(self.vertices[i][0]-self.center.x())**2)
            anglem=math.asin(anglem)
            
            if(self.vertices[i][0]>self.center.x() and self.vertices[i][2]>self.center.y()):
                self.angles.append(anglem)
            elif(self.vertices[i][0]>self.center.x() and self.vertices[i][2]<self.center.y()):
                self.angles.append(2*math.pi+anglem)
            elif (self.vertices[i][0]<self.center.x() and self.vertices[i][2]>self.center.y()):
                self.angles.append(math.pi-anglem)
            elif (self.vertices[i][0]<self.center.x() and self.vertices[i][2]<self.center.y()):
                self.angles.append(math.pi-anglem)
                
        
    def setLineWidth(self,width):
        self.lineWidth=width
        self.update()
        
        
    def paintEvent(self,evt):
        painter=QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        painter.setBrush(Canvas.backgroundColor)
        painter.drawRect(self.rect())
        
        painter.setBrush(QtCore.Qt.NoBrush)# no fill inside the shape

        if self.originalEllipse==True:
            self.draw_originalEllipse(painter)
        
        if self.standardEllipse==True:
            self.draw_standardEllipse(painter)
                    
        if self.generalisedEllipse==True:
            self.draw_generalisedEllipse(painter)
             
               
    def draw_generalisedEllipse(self,painter):
        penColor=QtGui.QColor(0,80,0) 
        pen=QtGui.QPen()
        pen.setColor(penColor)
        pen.setWidthF(self.lineWidth)
        painter.setPen(pen)
        I=self.numPt
        J=0
        
        if self.single_piece_mode==True:
            if self.manualJ_mode==True:
                J=self.manualJ
            
            elif self.autoJ_mode==True:
                J=self.findJ()
                self.autoJ=J
            
            if self.manualJ_mode==True or self.autoJ_mode==True:
                self.a,self.b=self.getCoefficients(J)
                self.generalisedEllipseVertices,self.Ea,self.Em=self.formGeneralizedEllipse(self.a,self.b)
     
                for i in range(I):
                    painter.drawLine(self.generalisedEllipseVertices[i][0],self.generalisedEllipseVertices[i][1],self.generalisedEllipseVertices[(i+1)%I][0],self.generalisedEllipseVertices[(i+1)%I][1])    
        
        elif self.segment_mode==True:        
            if self.manualJ_mode==True:
                J=self.manualJ
                a_first_half,b_first_half=self.coefficients_solver_for_first_half_of_segmented_ellipse(J)
                a_second_half,b_second_half=self.coefficients_solver_for_second_half_of_segmented_ellipse(J,a_first_half,b_first_half)
                
            first_half_vertices,first_half_Ea,first_half_Em=self.form_vertices_of_first_half_of_segmented_ellipse(a_first_half,b_first_half)    
            
        
    def form_vertices_of_first_half_of_segmented_ellipse(self,a,b):
        #CoefficientMatrix
        I=self.numPt/2
        first_half_vertices=[[0 for i in range(2)] for j in range(I)] 
        Ea=0.0
        Em=0.0
        d=[]
        J=len(a)/2
        for i in range(I):
            first_half_vertices[i][0]=self.center_first_half.x()+a[0]
            first_half_vertices[i][1]=self.center_first_half.y()+b[0]
            v=self.angles[i]
            for j in range(1,J+1):
                first_half_vertices[i][0]+=a[2*j-1]*math.cos(j*v)+a[2*j]*math.sin(j*v)
                first_half_vertices[i][1]+=b[2*j-1]*math.sin(j*v)+b[2*j]*math.cos(j*v)
           
            di=math.sqrt((self.vertices_first_half[i][0]-first_half_vertices[i][0])**2+(self.vertices_first_half[i][2]-first_half_vertices[i][1])**2)
            d.append(di)
            Ea+=(di/self.d_bar[i])
            if Em<di/self.d_bar[i]:
                Em=di/self.d_bar[i]
        Ea=Ea/(I+1)
        
        return generalisedEllipseVertices,Ea,Em
        
                    
    def coefficients_solver_for_first_half_of_segmented_ellipse(self,J):
        aConstArray=np.zeros(2*J+1)
        aCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')# row-major
        
        bConstArray=np.zeros(2*J+1)
        bCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')
        
        for i in range(I):
            aCoefficientMatrix[0,i]=1.
            bCoefficientMatrix[0,i]=1.
            
        for i in range(I):# for aCoefficientMatrix's column
            for j in range(1,J+1):# for aCoefficientMatrix's row
                try:
                    vi=self.angles_first_half[i]
                except IndexError:
                    error_dialog = QtWidgets.QErrorMessage(self)
                    error_dialog.showMessage('vertices data is empty, please choose a data file first')
                else:    
                    aCoefficientMatrix[2*j-1,i]=math.cos(vi*j)
                    aCoefficientMatrix[2*j,i]=math.sin(vi*j)

                    # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                    aConstArray[2*j-1]+=(self.vertices_first_half[i][0]-self.center_first_half.x())*math.cos(vi*j)
                    aConstArray[2*j]+=(self.vertices_first_half[i][0]-self.center_first_half.x())*math.sin(vi*j)
                
                    bCoefficientMatrix[2*j-1,i]=math.sin(vi*j)
                    bCoefficientMatrix[2*j,i]=math.cos(vi*j)
                   
                    bConstArray[2*j-1]+=(self.vertices_first_half[i][2]-self.center_first_half.y())*math.sin(vi*j)
                    bConstArray[2*j]+=(self.vertices_first_half[i][2]-self.center_first_half.y())*math.cos(vi*j)
                                  
        A=np.dot(aCoefficientMatrix,aCoefficientMatrix.transpose())
        a=np.linalg.solve(A,aConstArray)   
        B=np.dot(bCoefficientMatrix,bCoefficientMatrix.transpose())      
        b=np.linalg.solve(B,bConstArray) 
        
        return a,b
        
        
    def coefficients_solver_for_second_half_of_segmented_ellipse(self,J,a_first_half,b_first_half):
        aConstArray=np.zeros(2*J+1)
        aCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')# row-major
        
        bConstArray=np.zeros(2*J+1)
        bCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')
        
            
        for i in range(I):# for aCoefficientMatrix's column
            vi=self.angles_first_half[i]
           
            for j in range(3,J+1):# for aCoefficientMatrix's row
                aCoefficientMatrix[2*j-1,i]=math.cos(vi*j)
                aCoefficientMatrix[2*j,i]=math.sin(vi*j)
    
                # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                aConstArray[2*j-1]+=(self.vertices_first_half[i][0]-self.center_first_half.x())*math.cos(vi*j)
                aConstArray[2*j]+=(self.vertices_first_half[i][0]-self.center_first_half.x())*math.sin(vi*j)
                    
                        bCoefficientMatrix[2*j-1,i]=math.sin(vi*j)
                        bCoefficientMatrix[2*j,i]=math.cos(vi*j)
                       
                        bConstArray[2*j-1]+=(self.vertices_first_half[i][2]-self.center_first_half.y())*math.sin(vi*j)
                        bConstArray[2*j]+=(self.vertices_first_half[i][2]-self.center_first_half.y())*math.cos(vi*j)
                                      
        A=np.dot(aCoefficientMatrix,aCoefficientMatrix.transpose())
        a=np.linalg.solve(A,aConstArray)   
        B=np.dot(bCoefficientMatrix,bCoefficientMatrix.transpose())      
        b=np.linalg.solve(B,bConstArray) 
        
        return a,b
        
       
    def findJ(self):
        J=0
        a3,b3=self.getCoefficients(3)
        v3,Ea3,Em3=self.formGeneralizedEllipse(a3,b3)
        a10,b10=self.getCoefficients(10)
        v10,Ea10,Em10=self.formGeneralizedEllipse(a10,b10)
        if Ea3<Canvas.Ea_criteria and Em3<Canvas.Em_criteria:
            J=self.find_smaller_J(3,10,Ea3,Ea10)
        elif Ea10>=Canvas.Ea_criteria or Em10>=Canvas.Em_criteria:
            J=self.find_bigger_J(3,10,Ea3,Ea10,Canvas.Ea_criteria)
        elif Ea3>=Canvas.Ea_criteria or Em3>=Canvas.Em_criteria:
            J=self.find_inbetween_J(3,10,Ea3,Ea10,Em3,Em10)       
        return J
    
    
    def find_inbetween_J(self,J_small,J_big,Ea_smallJ,Ea_bigJ,Em_smallJ,Em_bigJ):
        #Linear interpolate to find J_small<J<J_big  
        
        if J_small>J_big:
            tmp=J_small
            J_small=J_big
            J_big=tmp
        
        if Ea_smallJ<Canvas.Ea_criteria and Ea_bigJ>=Canvas.Ea_criteria:
            tmp=Ea_smallJ
            Ea_smallJ=Ea_bigJ
            Ea_bigJ=tmp
            
        if Em_smallJ<Canvas.Em_criteria and Em_bigJ>=Canvas.Em_criteria:
            tmp=Em_smallJ
            Em_smallJ=Em_bigJ
            Em_bigJ=tmp  
        
        if J_small==J_big-1:
            return J_big
        
        J=0
        if Ea_smallJ>=Canvas.Ea_criteria:          
            J=int((Canvas.Ea_criteria-Ea_smallJ)*(J_big-J_small)/(Ea_bigJ-Ea_smallJ))+J_small
        elif Em_smallJ>=Canvas.Em_criteria:
            J=int((Canvas.Em_criteria-Em_smallJ)*(J_big-J_small)/(Em_bigJ-Em_smallJ))+J_small
        
        if J==J_small:
            J+=1
            
        a,b=self.getCoefficients(J)       
        v,Ea,Em=self.formGeneralizedEllipse(a,b)
        if Ea<Canvas.Ea_criteria and Em<Canvas.Em_criteria:
            if J==J_small+1:
                return J
            else:
                return self.find_inbetween_J(J_small,J,Ea_smallJ,Ea,Em_smallJ,Em)
        elif Ea>=Canvas.Ea_criteria or Em>=Canvas.Em_criteria:
            if J==J_big-1:
                return J_big
            else:
                return self.find_inbetween_J(J,J_big,Ea,Ea_bigJ,Em,Em_bigJ)
       
       
    def find_bigger_J(self,J_small,J_big,Ea_smallJ,Em_smallJ,Ea_bigJ,Em_bigJ):
        #Linear extrapolate to find bigger J>J_small
        if Ea_bigJ>=Canvas.Ea_criteria:
            J=int((Ea_criteria-Ea_smallJ)*(J_big-J_small)/(Ea_bigJ-Ea_smallJ))+J_small
        elif Em_bigJ>=Canvas.Em_criteria:
            J=int((Em_criteria-Em_smallJ)*(J_big-J_small)/(Em_bigJ-Em_smallJ))+J_small            
        a,b=self.getCoefficients(J)
        v,Ea,Em=self.formGeneralizedEllipse(a,b)
        if Ea>=Canvas.Ea_criteria or Em>=Canvas.Em_criteria: 
            return self.find_bigger_J(J,J_big,Ea,Ea_bigJ,Canvas.Ea_criteria)
        else:
            # we are close to the solution, hence, a while function will suffice
            while Ea<Canvas.Ea_criteria and Em<Canvas.Em_criteria and J>J_small:
                J-=1
                a,b=self.getCoefficients(J)
                v,Ea,Em=self.formGeneralizedEllipse(a,b)
            
            return J+1
    

    def find_smaller_J(self,J_small,J_big,Ea_smallJ,Ea_bigJ):
        #Linear extrapolate to find smaller J<J_small<J_big,
        #The criteria is always Canvas.Ea_criteria, 
        #because both (Ea_smallJ,Em_smallJ) and (Ea_bigJ,Em_bigJ) meet criteria. 
        #In this case, we will always use the average error Ea for the extrapolation since the average error is a global measurement. 
       
        J=int((Canvas.Ea_criteria-Ea_smallJ)*(J_big-J_small)/(Ea_bigJ-Ea_smallJ))+J_small
        a,b=self.getCoefficients(J)
        v,Ea,Em=self.formGeneralizedEllipse(a,b)
        if Ea<Canvas.Ea_criteria and Em<Canvas.Em_criteria:
            return self.fingSmallerJ(J,J_small,Ea,Ea_smallJ)
        else:
            # we are close to the solution, hence, a while function will suffice
            while Ea>=Canvas.Ea_criteria or Em>=Canvas.Em_criteria and J<J_small:
                J+=1
                a,b=self.getCoefficients(J)
            
            return J
                  
    
    def getCoefficients(self,J):#abtain a[2j+1] and b[2j+1]
        I=self.numPt
        aConstArray=np.zeros(2*J+1)
        aCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')# row-major
        
        bConstArray=np.zeros(2*J+1)
        bCoefficientMatrix=np.ndarray(shape=(2*J+1,I), dtype=float, order='C')
        
        for i in range(I):
            aCoefficientMatrix[0,i]=1.
            bCoefficientMatrix[0,i]=1.
            
        for i in range(I):# for aCoefficientMatrix's column
            for j in range(1,J+1):# for aCoefficientMatrix's row
                try:
                    vi=self.angles[i]
                except IndexError:
                    error_dialog = QtWidgets.QErrorMessage(self)
                    error_dialog.showMessage('Please choose a data file first')
                else:    
                    aCoefficientMatrix[2*j-1,i]=math.cos(vi*j)
                    aCoefficientMatrix[2*j,i]=math.sin(vi*j)

                    # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                    aConstArray[2*j-1]+=(self.vertices[i][0]-self.center.x())*math.cos(vi*j)
                    aConstArray[2*j]+=(self.vertices[i][0]-self.center.x())*math.sin(vi*j)
                
                    bCoefficientMatrix[2*j-1,i]=math.sin(vi*j)
                    bCoefficientMatrix[2*j,i]=math.cos(vi*j)
                   
                    bConstArray[2*j-1]+=(self.vertices[i][2]-self.center.y())*math.sin(vi*j)
                    bConstArray[2*j]+=(self.vertices[i][2]-self.center.y())*math.cos(vi*j)
                                  
        A=np.dot(aCoefficientMatrix,aCoefficientMatrix.transpose())
        a=np.linalg.solve(A,aConstArray)   
        B=np.dot(bCoefficientMatrix,bCoefficientMatrix.transpose())      
        b=np.linalg.solve(B,bConstArray) 
        
        return a,b
        
    
    def formGeneralizedEllipse(self,a,b):
        #CoefficientMatrix
        I=self.numPt
        generalisedEllipseVertices=[[0 for i in range(2)] for j in range(I)] 
        Ea=0.0
        Em=0.0
        d=[]
        J=len(a)/2
        for i in range(I):
            generalisedEllipseVertices[i][0]=self.center.x()+a[0]
            generalisedEllipseVertices[i][1]=self.center.y()+b[0]
            v=self.angles[i]
            for j in range(1,J+1):
                generalisedEllipseVertices[i][0]+=a[2*j-1]*math.cos(j*v)+a[2*j]*math.sin(j*v)
                generalisedEllipseVertices[i][1]+=b[2*j-1]*math.sin(j*v)+b[2*j]*math.cos(j*v)
           
            di=math.sqrt((self.vertices[i][0]-generalisedEllipseVertices[i][0])**2+(self.vertices[i][2]-generalisedEllipseVertices[i][1])**2)
            d.append(di)
            Ea+=(di/self.d_bar[i])
            if Em<di/self.d_bar[i]:
                Em=di/self.d_bar[i]
        Ea=Ea/self.numPt
        
        return generalisedEllipseVertices,Ea,Em
             
   
    def draw_originalEllipse(self,painter):
        penColor=QtCore.Qt.red   
        pen=QtGui.QPen()
        pen.setColor(penColor)
        pen.setWidthF(self.lineWidth)
        painter.setPen(pen)

        try:
            startPoint=QtCore.QPointF(self.vertices[0][0],self.vertices[0][2])
        except IndexError:
            error_dialog = QtWidgets.QErrorMessage(self)
            error_dialog.showMessage('Please choose a data file first')
        else:
            for i in range(self.numPt):
                p1=QtCore.QPointF(self.vertices[i-1][0],self.vertices[i-1][2])
                p2=QtCore.QPointF(self.vertices[i][0],self.vertices[i][2])
                painter.drawLine(p1,p2)
            """
            painterPath=QtGui.QPainterPath(startPoint)
            for i in range(self.numPt):
                tmp=self.vertices[(i+1)%self.numPt]-self.vertices[i-1]
                tmp/=math.sqrt(tmp[0]*tmp[0]+tmp[2]*tmp[2])#normalise vector
                arc=self.vertices[(i+1)%self.numPt]-self.vertices[i]
                arc=math.sqrt(arc[0]*arc[0]+arc[2]*arc[2])
                tmp=tmp*arc/1.
                controlPt1=tmp+self.vertices[i]
                    
                tmp=self.vertices[i]-self.vertices[(i+2)%self.numPt]
                tmp/=math.sqrt(tmp[0]*tmp[0]+tmp[2]*tmp[2])#normalise vector
                tmp=tmp*arc/1.
                controlPt2=tmp+self.vertices[(i+1)%self.numPt]
                painterPath.cubicTo(controlPt1[0],controlPt1[2],controlPt2[0],controlPt2[2],self.vertices[(i+1)%self.numPt][0],self.vertices[(i+1)%self.numPt][2])
                painterPath.moveTo(self.vertices[(i+1)%self.numPt][0],self.vertices[(i+1)%self.numPt][2])
                
            painter.drawPath(painterPath)    
            """          
                
        
    def draw_standardEllipse(self,painter): 
        penColor=QtCore.Qt.blue        
        pen=QtGui.QPen()
        pen.setColor(penColor)
        pen.setWidthF(self.lineWidth)
        painter.setPen(pen)
        
        # calculate the major axis and minor axis
        a0u=0.0
        a0b=0.0
        b0u=0.0
        b0b=0.0
        try:
            for i in range(self.numPt):
                angle=i*2*math.pi/self.numPt
                a0u=a0u+(self.vertices[i][0]-self.center.x())*math.cos(angle)
                a0b=a0b+math.cos(angle)*math.cos(angle)
                b0u=b0u+(self.vertices[i][2]-self.center.y())*math.sin(angle)
                b0b=b0b+math.sin(angle)*math.sin(angle)

        except IndexError:
            error_dialog = QtWidgets.QErrorMessage(self)
            error_dialog.showMessage('Please choose a data file first')
        else:
            width=a0u/a0b
            height=b0u/b0b
            
            painter.drawEllipse(self.center,width,height)


    def getEm(self):
        return self.Em
        
        
    def getEa(self):
        return self.Ea
        

if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("generalizedEllipseWin",exists=True):
        cmds.deleteUI("generalizedEllipseWin",wnd=True)
        
    w=CurveFittingWindowUI()
    w.show()  
    #w.raise_() 

