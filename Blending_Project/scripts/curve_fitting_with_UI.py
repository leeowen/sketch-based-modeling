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
        self.Ea_lineEdit.setMaxLength(5)

        self.Em_lineEdit=QtWidgets.QLineEdit()
        self.Em_lineEdit.setFixedWidth(200)
        self.Em_lineEdit.setMaxLength(5)
        
        self.canvas=Canvas()     
        
        self.lineWidth_lineEdit=QtWidgets.QLineEdit()
          
        
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
        
        self.J_spinBox.setVisible(False)
        self.Ea_lineEdit.setVisible(False)
        self.Em_lineEdit.setVisible(False)
        self.segment_comboBox.setVisible(False)
        self.saveAsDat_button.setVisible(False)
        self.saveAsImg_button.setVisible(False)
        

    def create_layout(self):   
        """     
        widget_layout_central=QtWidgets.QWidget()
        self.setCentralWidget(widget_layout_central)
        main_layout=QtWidgets.QHBoxLayout(widget_layout_central)
        """
        main_layout=QtWidgets.QHBoxLayout(self)
                 
        center_layout=QtWidgets.QHBoxLayout()
        center_layout.addWidget(self.canvas)              
        center_layout.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)        
              
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
        self.J_spinBox.valueChanged.connect(self.canvas.setJ)
        
            
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
        self.J_spinBox.setVisible(not checked)
        self.Ea_lineEdit.setVisible(not checked)
        self.Em_lineEdit.setVisible(not checked)
        self.segment_comboBox.setVisible(not checked)
        self.saveAsDat_button.setVisible(not checked)


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
        self.canvas.generalisedEllipse=checked
        self.canvas.update()
        s1=str(self.canvas.getEa()*100)
        s1=s1[0:4]
        s2=str(self.canvas.getEm()*100)
        s2=s2[0:4]
        self.Ea_lineEdit.setText(s1+'%')     
        self.Em_lineEdit.setText(s2+'%')
 
        
        self.J_spinBox.setVisible(checked)
        self.Ea_lineEdit.setVisible(checked)
        self.Em_lineEdit.setVisible(checked)
        self.segment_comboBox.setVisible(checked)
        self.saveAsDat_button.setVisible(checked)  
        self.saveAsImg_button.setVisible(checked)  
                  
        
    def drawMode_originalEllipse(self,checked):
        self.canvas.originalEllipse=checked
        self.canvas.update()

        self.saveAsDat_button.setVisible(checked)
        self.saveAsImg_button.setVisible(checked)
            
              
class Canvas(QtWidgets.QDialog):
    backgroundColor=QtCore.Qt.white
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
        self.J=10
        self.a=[]
        self.b=[]
        
    def sizeHint(self):
        return QtCore.QSize(650,600)
        
        
    def minimumSize(self):
        return QtCore.Qsize(650,600)
        
        
    def setJ(self,j):
        self.J=j
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
                
        
    def paintEvent(self,evt):
        painter=QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        painter.setBrush(Canvas.backgroundColor)
        painter.drawRect(self.rect())
        
        painter.setBrush(QtCore.Qt.NoBrush)# no fill inside the shape
        
        if self.generalisedEllipse==True:
            self.draw_generalisedEllipse(painter)

        if self.originalEllipse==True:
            self.draw_originalEllipse(painter)
        
        if self.standardEllipse==True:
            self.draw_standardEllipse(painter)
             
               
    def draw_generalisedEllipse(self,painter):
        penColor=QtGui.QColor(0,80,0) 
        painter.setPen(penColor)
        I=self.numPt
        aConstArray=np.zeros(2*self.J+1)
        aCoefficientMatrix=np.ndarray(shape=(2*self.J+1,I), dtype=float, order='C')# row-major
        aTrignometricMatrix=np.ndarray(shape=(I,2*self.J+1), dtype=float, order='C')
        bConstArray=np.zeros(2*self.J+1)
        bCoefficientMatrix=np.ndarray(shape=(2*self.J+1,I), dtype=float, order='C')
        bTrignometricMatrix=np.ndarray(shape=(I,2*self.J+1), dtype=float, order='C')
        for i in range(I):
            aCoefficientMatrix[0,i]=1.
            aTrignometricMatrix[i,0]=1.
            bCoefficientMatrix[0,i]=1.
            bTrignometricMatrix[i,0]=1.
            
        for i in range(I):# for aCoefficientMatrix's column, and trignomatricMatrix's row
            for j in range(1,self.J+1):# for aCoefficientMatrix's row, and trignomatricMatrix's column
                try:
                    vi=self.angles[i]
                except IndexError:
                    error_dialog = QtWidgets.QErrorMessage(self)
                    error_dialog.showMessage('Please choose a data file first')
                else:    
                    aCoefficientMatrix[2*j-1,i]=math.cos(vi*j)
                    aCoefficientMatrix[2*j,i]=math.sin(vi*j)
                    aTrignometricMatrix[i,2*j-1]=math.cos(vi*j)
                    aTrignometricMatrix[i,2*j]=math.sin(vi*j)
                    # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                    aConstArray[2*j-1]+=(self.vertices[i][0]-self.center.x())*math.cos(vi*j)
                    aConstArray[2*j]+=(self.vertices[i][0]-self.center.x())*math.sin(vi*j)
                
                    bCoefficientMatrix[2*j-1,i]=math.sin(vi*j)
                    bCoefficientMatrix[2*j,i]=math.cos(vi*j)
                    bTrignometricMatrix[i,2*j-1]=math.sin(vi*j)
                    bTrignometricMatrix[i,2*j]=math.cos(vi*j)
                    bConstArray[2*j-1]+=(self.vertices[i][2]-self.center.y())*math.sin(vi*j)
                    bConstArray[2*j]+=(self.vertices[i][2]-self.center.y())*math.cos(vi*j)
                                  
        A=np.dot(aCoefficientMatrix,aTrignometricMatrix)
        self.a=np.linalg.solve(A,aConstArray)   
        B=np.dot(bCoefficientMatrix,bTrignometricMatrix)      
        self.b=np.linalg.solve(B,bConstArray)           

        #CoefficientMatrix
        generalisedEllipseVertices=[[0 for i in range(2)] for j in range(I)] 
        self.Ea=0.0
        self.Em=0.0
        for i in range(I):
            generalisedEllipseVertices[i][0]=self.center.x()+self.a[0]
            generalisedEllipseVertices[i][1]=self.center.y()+self.b[0]
            v=self.angles[i]
            for j in range(1,self.J+1):
                generalisedEllipseVertices[i][0]+=self.a[2*j-1]*math.cos(j*v)+self.a[2*j]*math.sin(j*v)
                generalisedEllipseVertices[i][1]+=self.b[2*j-1]*math.sin(j*v)+self.b[2*j]*math.cos(j*v)
           
            di=math.sqrt((self.vertices[i][0]-generalisedEllipseVertices[i][0])**2+(self.vertices[i][2]-generalisedEllipseVertices[i][1])**2)
            self.d.append(di)
            self.Ea+=(di/self.d_bar[i])
            if self.Em<di/self.d_bar[i]:
                self.Em=di/self.d_bar[i]
        self.Ea=self.Ea/self.numPt
        
        for i in range(I):
            painter.drawLine(generalisedEllipseVertices[i][0],generalisedEllipseVertices[i][1],generalisedEllipseVertices[(i+1)%I][0],generalisedEllipseVertices[(i+1)%I][1])    
               
    
    def draw_originalEllipse(self,painter):
        penColor=QtCore.Qt.red   
        painter.setPen(penColor)

        try:
            startPoint=QtCore.QPointF(self.vertices[0][0],self.vertices[0][2])
        except IndexError:
            error_dialog = QtWidgets.QErrorMessage(self)
            error_dialog.showMessage('Please choose a data file first')
        else:
            """
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
                
        
    def draw_standardEllipse(self,painter): 
        penColor=QtCore.Qt.blue        
        painter.setPen(penColor)
        
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

