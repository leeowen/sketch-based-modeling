import maya.api.OpenMaya as om
import maya.cmds as cmds
from PySide2 import QtCore, QtGui, QtWidgets
from shiboken2 import wrapInstance
import maya.OpenMayaUI as omui
import math, sys, os, random
sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
sys.path.append(cmds.workspace(fn=True)+'/scripts/')
import numpy as np
import curve_fitting


def maya_main_window():
    main_window_ptr=omui.MQtUtil.mainWindow()
    return wrapInstance(long(main_window_ptr), QtWidgets.QWidget)
    
             
class CurveFittingWindowUI(QtWidgets.QWidget):
    
    def __init__(self, parent=maya_main_window()):
        super(CurveFittingWindowUI, self).__init__(parent)
        
        self.setWindowTitle("Generalized Ellipse")
        self.setMinimumSize(900, 700)
        self.setObjectName("generalizedEllipseWin")
        self.setWindowFlags(QtCore.Qt.WindowType.Dialog)
                       
        self.create_widgets()   
        self.create_layout()
        self.create_connection()
                                   
    
    def sizeHint(self):
        return QtCore.QSize(900, 700)
        
        
    def minimumSize(self):
        return QtCore.Qsize(900, 700)
        
        
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

        self.showPointIndex_checkBox = QtWidgets.QCheckBox('show point index')
        self.showPointIndex_checkBox.setChecked(False)
        self.showPointIndex_checkBox.setVisible(False)

        self.showTangent_checkBox = QtWidgets.QCheckBox('show tangent')
        self.showTangent_checkBox.setChecked(False)
        self.showTangent_checkBox.setVisible(False)

        self.showCurvature_checkBox = QtWidgets.QCheckBox('show curvature')
        self.showCurvature_checkBox.setChecked(False)
        self.showCurvature_checkBox.setVisible(False)

        self.showCutPointCandidate_checkBox = QtWidgets.QCheckBox('show cut point candidate(s)')
        self.showCutPointCandidate_checkBox.setChecked(False)
        self.showCutPointCandidate_checkBox.setVisible(False)

        self.cut_point_label = QtWidgets.QLabel('cut point(s) list:')
        self.cut_point_lineEdit = QtWidgets.QLineEdit()
        self.cut_point_button = QtWidgets.QPushButton('Confirm')
                       
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
        self.J_spinBox.setMinimum(3)
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
        self.segment_comboBox.addItem('2 symmetrical halves')
        self.segment_comboBox.addItem('fragment')
        self.segment_comboBox.addItem('composite')
        
        self.range_label=QtWidgets.QLabel('range:')
        self.range_label.setMaximumWidth(150)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        self.range_label.setSizePolicy(sizePolicy)
        self.range_slider=QtWidgets.QSlider(QtCore.Qt.Horizontal)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.range_slider.setSizePolicy(sizePolicy)
        self.range_slider.setMaximumWidth(150)
        self.range_slider.setMinimum(0)
        self.range_slider.setMaximum(100)
        
        self.radio_group.setVisible(False)
        self.segment_label.setVisible(False)
        self.segment_comboBox.setVisible(False)
        self.range_label.setVisible(False)
        self.range_slider.setVisible(False)
        self.saveAsDat_button.setVisible(False)
        self.saveAsImg_button.setVisible(False)
        self.cut_point_button.setVisible(False)
        self.cut_point_label.setVisible(False)
        self.cut_point_lineEdit.setVisible(False)
        

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

        checkBox2_layout=QtWidgets.QVBoxLayout()
        checkBox2_layout.addWidget(self.showPointIndex_checkBox)
        checkBox2_layout.addWidget(self.showTangent_checkBox)
        checkBox2_layout.addWidget(self.showCurvature_checkBox)
        checkBox2_layout.addWidget(self.showCutPointCandidate_checkBox)
        right_layout.addLayout(checkBox2_layout)

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
        
        rangeSlider_layout=QtWidgets.QHBoxLayout()
        rangeSlider_layout.addWidget(self.range_label)
        rangeSlider_layout.addWidget(self.range_slider)
        right_layout.addLayout(rangeSlider_layout) 
        
        right_layout.addStretch(1)

        cutPoint_layout = QtWidgets.QVBoxLayout()
        cutPoint_layout.addWidget(self.cut_point_label)
        cutPoint_layout.addWidget(self.cut_point_lineEdit)
        cutPoint_layout.addWidget(self.cut_point_button)
        right_layout.addLayout(cutPoint_layout)

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
        self.range_slider.valueChanged.connect(self.fragment_range_change)
        self.showPointIndex_checkBox.toggled.connect(self.update_display_point_index)
        self.showCurvature_checkBox.toggled.connect(self.update_curvature)
        self.showTangent_checkBox.toggled.connect(self.update_tangent)
        self.cut_point_button.clicked.connect(self.get_cut_point)
        
    
    def fragment_range_change(self,value):
        self.canvas.set_fragment_range(value/100.0)
        
        
    def segmentMode_change(self,text):
        if text=='single piece':
            self.canvas.single_piece_mode=True
            self.canvas.symmetry_mode=False
            self.canvas.fragment_mode=False
            self.canvas.composite_mode=False
            self.range_label.setVisible(False)
            self.range_slider.setVisible(False)
            self.manualJ_mode_radioButton.setVisible(True)
            self.cut_point_button.setVisible(False)
            self.cut_point_label.setVisible(False)
            self.cut_point_lineEdit.setVisible(False)
        elif text=='2 symmetrical halves':
            self.canvas.symmetry_mode=True
            self.canvas.single_piece_mode=False
            self.canvas.fragment_mode=False
            self.canvas.composite_mode=False
            self.range_label.setVisible(False)
            self.range_slider.setVisible(False)
            self.manualJ_mode_radioButton.setVisible(True)
            self.cut_point_button.setVisible(False)
            self.cut_point_label.setVisible(False)
            self.cut_point_lineEdit.setVisible(False)
        elif text=='fragment':
            self.canvas.symmetry_mode = False
            self.canvas.single_piece_mode = False
            self.canvas.fragment_mode = True
            self.canvas.composite_mode = False
            self.autoJ_mode_radioButton.setVisible(False)
            self.range_label.setVisible(True)
            self.range_slider.setVisible(True)
            self.manualJ_mode_radioButton.setVisible(True)
            self.cut_point_button.setVisible(False)
            self.cut_point_label.setVisible(False)
            self.cut_point_lineEdit.setVisible(False)
        elif text == 'composite':
            self.canvas.symmetry_mode=False
            self.canvas.single_piece_mode=False
            self.canvas.fragment_mode=False
            self.canvas.composite_mode=True
            self.range_label.setVisible(False)
            self.range_slider.setVisible(False)
            self.manualJ_mode_radioButton.setVisible(True)
            self.autoJ_mode_radioButton.setVisible(True)
            self.showTangent_checkBox.setVisible(True)
            self.cut_point_button.setVisible(True)
            self.cut_point_label.setVisible(True)
            self.cut_point_lineEdit.setVisible(True)

            
        self.canvas.update()
        self.showEaEm()
                        
    
    def update_visibility_manual_J_mode(self,checked):
        self.canvas.manualJ_mode=checked
        self.canvas.autoJ_mode=not checked
        self.J_spinBox.setReadOnly(not checked)
        self.canvas.repaint()# update() prevents multiple fast repaints. call repaint() instead.
        self.showEaEm()
        if self.canvas.manualJ_value > self.J_spinBox.value():
            self.J_spinBox.setValue(self.canvas.manualJ_value)
        
        
    def update_visibility_auto_J_mode(self,checked):
        self.canvas.autoJ_mode=checked
        self.canvas.manualJ_mode=not checked
        self.J_spinBox.setReadOnly(checked)
        self.canvas.repaint()
        self.showEaEm()
        self.J_spinBox.setValue(self.canvas.autoJ_value)
        
    
    def update_manual_J_value(self,J):
        self.canvas.setManualJ(J)
        self.canvas.repaint()
        self.showEaEm()


    def update_display_point_index(self,checked):
        self.canvas.show_point_index = checked
        self.canvas.repaint()
        

    def update_tangent(self,checked):
        self.canvas.activateTangent = checked
        self.canvas.repaint()


    def update_curvature(self,checked):
        self.canvas.activateCurvature=checked
        self.canvas.repaint()


    def save_dataFn(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        dirPath+='data/'
        source_file_name = self.filePath_lineEdit.text()
        file_name = source_file_name.split('/')[-1]
        file_name = file_name.replace('Source_','')
        dirPath += file_name.split('_cross_section_')[0]
        directory = os.path.dirname(dirPath+'/'+ file_name)
        
        try:
            os.stat(directory)
        except:
            os.mkdir(directory) 
            
        FILE_FILTERS="DAT(*.dat);All Files(*.*)"
        selected_filter="DAT(*.dat)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_path,selected_filter=QtWidgets.QFileDialog.getSaveFileName(self, 'save',dirPath +'/'+ file_name,FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            f=open(file_path,"w+")
            f.write('range:')
            if self.canvas.single_piece_mode==True:
                f.write('0-360 \n')
            f.write('center: {} {}\n'.format(self.canvas.center[0],self.canvas.center[2]))
            f.write('a: ')
            for i in self.canvas.a:
                f.write(str(i)+' ')
            f.write('\n')
            f.write('b: ')
            for i in self.canvas.b:
                f.write(str(i)+' ')
            f.write('\n')
            #f.write('angles: ')
            #for i in self.canvas.angles:
            #    f.write(str(i)+' ')
            #f.write('\n')
            f.close()
        
        
    def save_imageFn(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        dirPath+='images/'
        FILE_FILTERS="PNG(*.png);All Files(*.*)"
        selected_filter="PNG(*.png)"# default filter, also store last selected filter and can be used as the default filter for next select
        image_name = self.filePath_lineEdit.text()
        if 'Source_' in image_name:
            image_name = image_name.split('Source_')[1]
            image_name = image_name.replace('.dat','.png')
            part_name = image_name.split('_cross_section')[0]
            image_name = image_name.replace('_cross_section','_generalised_ellipse')
            dirPath += part_name+ '/'
        try:
            os.stat(dirPath)
        except:
            os.mkdir(dirPath)
        file_path,selected_filter = QtWidgets.QFileDialog.getSaveFileName(self, 'save', dirPath+image_name, FILE_FILTERS, selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_path:
            a=file_path.find('.png')
            b=file_path.find('.PNG')
            if a==-1 and b==-1:
                file_path+='.png'
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
        self.save_file_and_image_toggle()


    def save_file_and_image_toggle(self):
        if self.originalEllipse_checkBox.isChecked() or self.standardEllipse_checkBox.isChecked() or self.generalizedEllipse_checkBox.isChecked():
            self.saveAsDat_button.setVisible(True)
            self.saveAsImg_button.setVisible(True)
        else:
            self.saveAsDat_button.setVisible(False)
            self.saveAsImg_button.setVisible(False)


    def update_visibility_standardEllipse_mode(self,checked):
        self.save_file_and_image_toggle()
        
    
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
        self.showPointIndex_checkBox.setVisible(checked)
        self.showTangent_checkBox.setVisible(checked)
        self.radio_group.setVisible(checked)
        self.segment_label.setVisible(checked)
        self.segment_comboBox.setVisible(checked)
        self.canvas.generalisedEllipse=checked
        self.save_file_and_image_toggle()
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

        self.showPointIndex_checkBox.setVisible(checked)
        self.showTangent_checkBox.setVisible(checked)
        self.save_file_and_image_toggle()
        self.showCurvature_checkBox.setVisible(checked)
        self.showTangent_checkBox.setVisible(checked)
        self.showCutPointCandidate_checkBox.setVisible(checked)


    def get_cut_point(self):
        str = self.cut_point_lineEdit.text()
        list = str.split(' ')
        cut_points = [int(i) for i in list]
        self.canvas.cut_points = cut_points
        self.canvas.vertices_matrix,self.canvas.angles_matrix,self.canvas.segment_center_list = curve_fitting.cut_curve(self.canvas.vertices, self.canvas.angles, self.canvas.cut_points, self.canvas.isClosed)
            
              
class Canvas(QtWidgets.QDialog):
    backgroundColor = QtCore.Qt.white
    Ea_criteria = 0.01
    Em_criteria = 0.025
    color_list = [QtCore.Qt.darkCyan, QtCore.Qt.magenta, QtCore.Qt.darkYellow, QtCore.Qt.darkRed,
                  QtCore.Qt.black, QtCore.Qt.darkMagenta]
    
    def __init__(self,parent=None):
        super(Canvas,self).__init__(parent)
        self.setMinimumSize(650,600)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.originalEllipse=False
        self.generalisedEllipse=False
        self.standardEllipse=False
        self.manualJ_value=20
        self.autoJ_value=0
        self.lineWidth=1
        self.manualJ_mode=False
        self.autoJ_mode=False
        self.single_piece_mode=True
        self.symmetry_mode=False
        self.fragment_mode=False
        self.composite_mode=False
        self.activateTangent=False
        self.activateCurvature=False
        self.activateCutPoint=False
        self.show_point_index=False
        self.isClosed=True

        
    def sizeHint(self):
        return QtCore.QSize(650,600)
        
        
    def minimumSize(self):
        return QtCore.Qsize(650,600)
        
        
    def setManualJ(self,j):
        self.manualJ_value=j
        self.update()


    def readFile(self,file_path):
        self.vertices = []
        self.angles = []
        self.d = []
        self.d_bar = []
        self.Ea = 0.0
        self.Em = 0.0
        self.numPt = 0
        self.center = om.MVector(0,0,0)
        self.fragment_range=0.0
        self.start_index=0
        self.end_index=0
        self.a=[]
        self.b=[]
        self.generalisedEllipseVertices=[]
        self.cut_points = []
        self.composite_vertices = []
        self.composite_a = []
        self.composite_b = []
        self.input_file_path = file_path
        self.vertices_matrix = []
        self.angles_matrix = []
        self.segment_center_list = []
        f = open(file_path,'r')
        content = f.readlines()
        for line in content:
            p=line.split()
            self.vertices.append(om.MVector(float(p[0]), float(p[1]), float(p[2])))
        f.close()
        self.center = curve_fitting.getCenter_2D(self.vertices)

        def sort():
            if self.vertices[0][0] > self.center[0] and self.vertices[-1][0] < self.center[0]:
                return
            elif self.vertices[0][0] < self.center[0]:
                i = 0
                while self.vertices[i][0] < self.center[0]:
                    i += 1
                for d in range(i):
                    b = self.vertices.pop(0)
                    self.vertices.append(b)
                return
            elif self.vertices[-1][0] > self.center[0]:
                i = -1
                while self.vertices[i][0] > self.center[0]:
                    i += 1
                for d in range(i):
                    self.vertices.insert(self.numPt, self.vertices[0])
                    self.vertices.remove(self.vertices[0])
                return

        def modify_first_and_middle_vertices():
            N = self.vertices[-1]
            O = self.vertices[0]
            C = self.center
            if N[0] != O[0]:
                k = (O[2] - N[2]) / (O[0] - N[0])
                b = N[2] - k * N[0]
                new_z = k * C[0] + b
                k = (O[1] - N[1]) / (O[0] - N[0])
                b = N[1] - k * N[0]
                new_y = k * C[0] + b
                first_vertex = om.MVector(C[0], new_y, new_z)
                self.vertices.insert(0, first_vertex)
                self.vertices.pop(len(self.vertices)-1)
            else:
                raise ValueError('the first and last vertex x position is the same')

            e = 1
            for i in range(1, len(self.vertices)):
                if self.vertices[i][0] > C[0] and self.vertices[i+1][0] < C[0]:
                    e = i
                    break
            E = self.vertices[e]
            F = self.vertices[e+1]

            if e + 1 != len(self.vertices)/2:
                cmds.confirmDialog(message="bad photo, center is not in the middle, don't use symmetry function", button=["ok"])

            if E[0] != F[0]:
                k = (E[2] - F[2]) / (E[0] - F[0])
                b = F[2] - k * F[0]
                new_z = k * C[0] + b
                k = (E[1] - F[1]) / (E[0] - F[0])
                b = E[1] - k * E[0]
                new_y = k * C[0] + b
                middle_vertex = om.MVector(C[0], new_y, new_z)
                self.vertices.insert(e+1, middle_vertex)
                self.vertices.pop(e+2)
            else:
                raise ValueError('two vertices in the middle have same x position')
            return e

        #sort()
        #I = modify_first_and_middle_vertices()
        self.numPt = len(self.vertices)
        self.d_bar = curve_fitting.get_d_bar_2D(self.vertices, self.center)
        self.angles = curve_fitting.calculateAngle_2D(self.vertices, self.center)
        self.tangents = curve_fitting.calculateTangent_2D(self.vertices, self.angles)
        self.normals = curve_fitting.calculateNormal_2D(self.tangents)
        self.curvatures = curve_fitting.calculateCurvature_2D(self.tangents, self.angles)
        self.isClosed = True
        angle_difference = abs(self.angles[-1] - self.angles[0] - 2 * math.pi)
        if angle_difference > math.pi * 2.0:
            angle_difference -= math.pi * 2.0
        if angle_difference > 3 * math.pi / self.numPt:
            self.isClosed = False
        
        # split data for 2 symmetrical segments respectively
        self.a_first_half = []
        self.b_first_half = []
        I = len(self.vertices)/2
        self.vertices_first_half = self.vertices[0:I+1]
        self.angles_first_half = self.angles[0:I+1]
        self.vertices_second_half = self.vertices[I:self.numPt]
        self.vertices_second_half.append(self.vertices[0])
        self.angles_second_half = self.angles[I:self.numPt]
        self.angles_second_half.append(self.angles[0])
        self.center_first_half = curve_fitting.getCenter_2D(self.vertices_first_half)
        self.center_second_half = curve_fitting.getCenter_2D(self.vertices_second_half)

        if 'Chest' in self.input_file_path:# rebuild the curve
            # delete the unnecessary points
            delete_points_list = {'Source_Chest_cross_section_u_at_22_percentage.dat': [34, 37, 83, 86],
                                  'Source_Chest_cross_section_u_at_66_percentage.dat': [27, 39, 81, 93],
                                  'Source_Chest_cross_section_u_at_86_percentage.dat': [27, 42, 78, 93],
                                  'Source_Chest_cross_section_u_at_96_percentage.dat': [27, 42, 78, 93]}
            self.cut_points = [27, 42, 78, 93]
            #self.vertices = curve_fitting.rebuild_curve(file_path, self.vertices, delete_points_list)
            # delete the unnecessary points
            filepath = self.input_file_path.split('/')[-1]
            if filepath != 'Source_Chest_cross_section_u_at_0_percentage_worldspace.dat':
                del self.vertices[delete_points_list.get(filepath)[2] + 1:delete_points_list.get(filepath)[3]]
                del self.vertices[delete_points_list.get(filepath)[0] + 1:delete_points_list.get(filepath)[1]]
                self.cut_points[1] -= delete_points_list.get(filepath)[1] - delete_points_list.get(filepath)[0] - 1
                self.cut_points[2] -= delete_points_list.get(filepath)[1] - delete_points_list.get(filepath)[0] - 1
                self.cut_points[3] -= delete_points_list.get(filepath)[1] - delete_points_list.get(filepath)[0] - 1 + \
                                 delete_points_list.get(filepath)[3] - delete_points_list.get(filepath)[2] - 1

            self.center = curve_fitting.getCenter_3D(self.vertices)
            self.angles = curve_fitting.calculateAngle_3D(self.vertices, self.center)
            self.d_bar = curve_fitting.get_d_bar_3D(self.vertices, self.center)
            self.numPt = len(self.vertices)

            # cut curve
            self.vertices_matrix, self.angles_matrix, self.segment_center_list = curve_fitting.cut_curve(self.vertices, self.angles, self.cut_points, self.isClosed)


    def setLineWidth(self,width):
        self.lineWidth = width
        self.update()
        
    
    def set_fragment_range(self,value):   
        self.fragment_range = value
        self.update()
        
         
    def paintEvent(self, evt):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
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

        if self.show_point_index==True:
            self.draw_point_index(painter)
            
        painter.end()
             

    def draw_point_index(self,painter):
        font = QtGui.QFont()
        font.setPointSize(5)
        painter.setFont(font)
        for i in range(self.numPt):
            p = self.vertices[i]
            painter.drawText(p[0] * 300 + self.width() / 2., p[2] * 300 + self.height() / 2., str(i))
            painter.drawEllipse(QtCore.QPointF(p[0] * 300 + self.width() / 2., p[2] * 300 + self.height() / 2.),1,1)


    def draw_generalisedEllipse(self, painter):
        penColor = QtGui.QColor(0, 80, 0)
        pen = QtGui.QPen()
        pen.setColor(penColor)
        pen.setWidthF(self.lineWidth)
        painter.setPen(pen)
        I = self.numPt
        J = 0

        if self.single_piece_mode == True:
            if self.manualJ_mode == True:
                J = self.manualJ_value
            
            elif self.autoJ_mode==True:
                J = curve_fitting.findJ_2D(self.vertices, self.angles, self.d_bar, self.center, self.Ea_criteria, self.Em_criteria, curve_fitting.getCoefficients_2D, curve_fitting.formGeneralizedEllipse_2D,0)
                self.autoJ_value = J

            if J != 0:
                self.a, self.b = curve_fitting.getCoefficients_2D(J, self.vertices, self.center, self.angles)
                self.generalisedEllipseVertices, self.Ea, self.Em = curve_fitting.formGeneralizedEllipse_2D(self.a, self.b, self.vertices, self.center, self.angles, self.d_bar)

                for i in range(I):
                    painter.drawLine(self.generalisedEllipseVertices[i][0]*300+self.width()/2., self.generalisedEllipseVertices[i][2]*300+self.height()/2.,
                                     self.generalisedEllipseVertices[(i+1)%I][0]*300+self.width()/2., self.generalisedEllipseVertices[(i+1) % I][2]*300+self.height()/2.)

            if self.activateTangent == True:
                pen1 = QtGui.QPen()
                pen1.setColor(QtCore.Qt.blue)
                pen1.setWidthF(self.lineWidth)
                painter.setPen(pen1)
                for i in range(I):
                    x_tan, y_tan, x, y = curve_fitting.position_and_tangent_of_parametric_point_2D(self.a, self.b, self.angles[i])
                    x += self.center[0]
                    y += self.center[2]
                    t = QtGui.QVector2D(x_tan, y_tan).normalized() * 20.0
                    p0 = QtCore.QPointF(x * 300 + self.width()/2., y * 300 + self.height()/2.)
                    p1 = QtCore.QPointF(x * 300 + self.width()/2. + t.x(), y * 300 + self.height()/2. + t.y())
                    painter.drawLine(p0, p1)
                painter.setPen(pen)

        elif self.symmetry_mode == True:
            if self.manualJ_mode == True:
                J = self.manualJ_value

            elif self.autoJ_mode==True:
                J = curve_fitting.findJ_2D(self.vertices_first_half, self.angles_first_half,self.d_bar,self.center_first_half,self.Ea_criteria,self.Em_criteria,
                                        curve_fitting.getCoefficients_2D,curve_fitting.form_vertices_of_fragment_2D,0)
                self.autoJ_value=J

            if self.manualJ_mode==True or self.autoJ_mode==True:
                a_first_half, b_first_half = curve_fitting.getCoefficients_2D(
                    J,self.vertices_first_half,self.center_first_half,self.angles_first_half)
                a_second_half, b_second_half = curve_fitting.getCoefficients_for_second_half_of_symmetrical_ellipse(
                    a_first_half, b_first_half)

                symmetry_ellipse_vertices, Ea, Em = curve_fitting.form_vertices_of_symmetry_ellipse_2D(
                    self.vertices_first_half,
                    self.center_first_half,
                    self.angles_first_half,
                    self.d_bar, a_first_half,
                    b_first_half, a_second_half,
                    b_second_half)
                self.Ea = Ea
                self.Em = Em
                for i in range(self.numPt / 2):
                    painter.drawLine(symmetry_ellipse_vertices[i][0] * 300 + self.width() / 2.,
                                     symmetry_ellipse_vertices[i][2] * 300 + self.height() / 2.,
                                     symmetry_ellipse_vertices[(i + 1) % self.numPt][0] * 300 + self.width() / 2.,
                                     symmetry_ellipse_vertices[(i + 1) % self.numPt][2] * 300 + self.height() / 2.)
                second_half_color = QtGui.QColor(127, 0, 255)
                pen2 = QtGui.QPen()
                pen2.setColor(second_half_color)
                pen2.setWidthF(self.lineWidth)
                painter.setPen(pen2)
                for i in range(self.numPt / 2, self.numPt):
                    painter.drawLine(symmetry_ellipse_vertices[i][0] * 300 + self.width() / 2.,
                                     symmetry_ellipse_vertices[i][2] * 300 + self.height() / 2.,
                                     symmetry_ellipse_vertices[(i + 1) % self.numPt][0] * 300 + self.width() / 2.,
                                     symmetry_ellipse_vertices[(i + 1) % self.numPt][2] * 300 + self.height() / 2.)
                painter.setPen(pen)

        elif self.fragment_mode==True:
            if self.manualJ_mode==True:
                J = self.manualJ_value
                angles, vertices, center, self.start_index, self.end_index = curve_fitting.extract_fragment_data(self.vertices, self.angles, self.fragment_range)
                a, b = curve_fitting.getCoefficients_for_fragmented_ellipse(J,angles,vertices,center)
                fragment_vertices, Ea, Em = curve_fitting.form_vertices_of_fragment_2D(a, b, vertices, center, angles, self.d_bar, self.start_index)
                self.Ea = Ea
                self.Em = Em
                for i in range(len(fragment_vertices)-1):
                    painter.drawLine(fragment_vertices[i][0]*300+self.width()/2., fragment_vertices[i][2]*300+self.height()/2., fragment_vertices[i+1][0]*300+self.width()/2., fragment_vertices[i+1][2]*300+self.height()/2.)

        elif self.composite_mode == True:
            J1 = 0
            self.J_total = 0

            if self.cut_points:
                # for the first segment
                self.composite_vertices = []
                self.composite_a = []
                self.composite_b = []

                J0 = curve_fitting.findJ_2D(self.vertices_matrix[0], self.angles_matrix[0], self.d_bar,
                                         self.segment_center_list[0], self.Ea_criteria, self.Em_criteria,curve_fitting.getCoefficients_2D,curve_fitting.form_vertices_of_fragment_2D,self.cut_points[0])

                self.J_total += J0
                print "J0 = {}".format(J0)
                a, b = curve_fitting.getCoefficients_2D(J0,self.vertices_matrix[0],self.segment_center_list[0],self.angles_matrix[0])

                vertices, self.Ea, self.Em = curve_fitting.form_vertices_of_fragment_2D(a, b, self.vertices_matrix[0],
                                                                                     self.segment_center_list[0],
                                                                                     self.angles_matrix[0], self.d_bar,
                                                                                     self.cut_points[0])
                self.composite_vertices.append(vertices)
                self.composite_a.append(a)
                self.composite_b.append(b)
            else:# do nothing because cut points haven't been provided
                return

            # in-between segments
            def inbetween_2D_open_curve_auto_mode():
                for i in range(len(self.cut_points)):
                    J, vertices, a, b, Ea, Em = curve_fitting.compisite_segment_with_one_end_shared_2D(i + 1,self.vertices_matrix,self.angles_matrix,self.composite_a,self.composite_b,self.segment_center_list,
                                      self.d_bar,self.Ea_criteria, self.Em_criteria,self.cut_points,self.isClosed)
                    cut_pt_index = self.cut_points[i]
                    self.Ea = (self.Ea * cut_pt_index + Ea * len(vertices)) / (cut_pt_index + len(vertices))
                    self.Em = max(self.Em, Em)
                    self.composite_vertices.append(vertices)
                    self.composite_a.append(a)
                    self.composite_b.append(b)
                    self.J_total += J
                    print "J{} = {}".format(i + 1, J)

            def inbetween_2D_closed_curve_auto_mode():
                for i in range(1, len(self.cut_points) - 1):
                    cut_pt_index = self.cut_points[i - 1]
                    J, vertices, a, b, Ea, Em = curve_fitting.compisite_segment_with_one_end_shared_2D(i, self.vertices_matrix,
                                                                                                       self.angles_matrix,
                                                                                                       self.composite_a,
                                                                                                       self.composite_b,
                                                                                                       self.segment_center_list,
                                                                                                       self.d_bar,
                                                                                                       self.Ea_criteria,
                                                                                                       self.Em_criteria,
                                                                                                       self.cut_points,
                                                                                                       self.isClosed)
                    self.Ea = (self.Ea * cut_pt_index + Ea * len(vertices)) / (cut_pt_index + len(vertices))
                    self.Em = max(self.Em, Em)
                    self.composite_vertices.append(vertices)
                    self.composite_a.append(a)
                    self.composite_b.append(b)
                    self.J_total += J

                    print "J{} = {}".format(i, J)

            if self.isClosed == False:  # curve is open
                inbetween_2D_open_curve_auto_mode()
            else:
                inbetween_2D_closed_curve_auto_mode()

            # for the end segment that links the first segment
            if self.isClosed == False:# because it is an open curve, there is no end segment links the first segment
                self.autoJ_value = self.J_total
                self.manualJ_value = self.J_total
            else:
                x0_tan, y0_tan, x0, y0 = curve_fitting.position_and_tangent_of_parametric_point_2D(
                    self.composite_a[-1],
                    self.composite_b[-1],
                    self.angles_matrix[-2][
                        -1])
                x0 += self.segment_center_list[-2][0]
                y0 += self.segment_center_list[-2][2]
                previous = {'position x': x0, 'position y': y0, 'tangent x': x0_tan, 'tangent y': y0_tan,
                            'cut point index': self.cut_points[-1]}
                x1_tan, y1_tan, x1, y1 = curve_fitting.position_and_tangent_of_parametric_point_2D(self.composite_a[0],
                                                                                                self.composite_b[0],
                                                                                                self.angles_matrix[0][0])
                x1 += self.segment_center_list[0][0]
                y1 += self.segment_center_list[0][2]
                next = {'position x': x1, 'position y': y1, 'tangent x': x1_tan, 'tangent y': y1_tan,
                        'cut point index': self.cut_points[0]}

                Jn = 0
                if self.manualJ_mode == True:
                    if self.manualJ_value - self.J_total < 3:
                        Jn = 3
                    else:
                        Jn = self.manualJ_value - self.J_total
                    self.J_total += Jn
                    self.manualJ_value = self.J_total
                elif self.autoJ_mode == True:
                    Jn = curve_fitting.findJ_for_end_segment_2D(self.vertices_matrix[-1],self.angles_matrix[-1], self.d_bar, self.segment_center_list[-1], self.Ea_criteria, self.Em_criteria, previous, next)
                    self.J_total += Jn
                    self.autoJ_value = self.J_total
                a, b = curve_fitting.getCoefficients_for_end_composite_2D(Jn, self.vertices_matrix[-1],
                                                                       self.segment_center_list[-1],
                                                                       self.angles_matrix[-1], previous, next)
                composite_vertices_n, Ean, Emn = curve_fitting.form_vertices_of_fragment_2D(a, b,
                                                                                         self.vertices_matrix[-1],
                                                                                         self.segment_center_list[-1],
                                                                                         self.angles_matrix[-1],
                                                                                         self.d_bar,
                                                                                         self.cut_points[-1])
                self.composite_a.append(a)
                self.composite_b.append(b)
                self.composite_vertices.append(composite_vertices_n)

                self.Ea = (self.Ea * (self.numPt - len(self.vertices_matrix[-1])) + Ean * len(
                    self.vertices_matrix[-1])) / self.numPt
                if self.Em < Emn:
                    self.Em = Emn
                self.manualJ_value = self.J_total + Jn
                print "J{} = {}".format(len(self.cut_points) - 1, Jn)
                def draw_meet_points_tangent():
                    pen1 = QtGui.QPen()
                    pen1.setColor(QtCore.Qt.magenta)
                    painter.setPen(pen1)
                    painter.drawLine(x0 * 300 + self.width() / 2., y0 * 300 + self.height() / 2., x0 * 300 + x0_tan * 10 + self.width() / 2., y0 * 300 + y0_tan * 10 + self.height() / 2.)
                    painter.drawLine(x1 * 300 + self.width() / 2., y1 * 300 + self.height() / 2., x1 * 300 + x1_tan * 10 + self.width() / 2., y1 * 300 + y1_tan * 10 + self.height() / 2.)
                    painter.setPen(pen)

                #draw_meet_points_tangent()

                #test_x0_tan, test_y0_tan, test_x0, test_y0 = curve_fitting.position_and_tangent_of_parametric_point_2D(self.composite_a[1], self.composite_b[1], self.angles_matrix[1][0])
                #test_x1_tan, test_y1_tan, test_x1, test_y1 = curve_fitting.position_and_tangent_of_parametric_point_2D(self.composite_a[1], self.composite_b[1], self.angles_matrix[1][-1])
                #test_x0 += self.segment_center_list[1][0]
                #test_y0 += self.segment_center_list[1][2]
                #test_x1 += self.segment_center_list[1][0]
                #test_y1 += self.segment_center_list[1][2]

                def test_meet_points_tangent():
                    pen1 = QtGui.QPen()
                    pen1.setColor(QtCore.Qt.black)
                    painter.setPen(pen1)
                    painter.drawLine(test_x0 * 300 + self.width() / 2., test_y0 * 300 + self.height() / 2., test_x0 * 300 + test_x0_tan * 10 + self.width() / 2., test_y0 * 300 + test_y0_tan * 10 + self.height() / 2.)
                    painter.drawLine(test_x1 * 300 + self.width() / 2., test_y1 * 300 + self.height() / 2., test_x1 * 300 + test_x1_tan * 10 + self.width() / 2., test_y1 * 300 + test_y1_tan * 10 + self.height() / 2.)
                    painter.setPen(pen)

                #test_meet_points_tangent()
            for i in range(0, len(self.composite_vertices)):
                penx = QtGui.QPen()
                penx.setColor(self.color_list[i])
                penx.setWidthF(self.lineWidth)
                painter.setPen(penx)
                row = self.composite_vertices[i]
                for j in xrange(len(row)-1):
                    painter.drawLine(row[j][0] * 300 + self.width() / 2., row[j][2] * 300 + self.height() / 2., row[j + 1][0] * 300 + self.width() / 2., row[j + 1][2] * 300 + self.height() / 2.)
            painter.setPen(pen)


    def draw_originalEllipse(self,painter):
        penColor=QtCore.Qt.red   
        pen=QtGui.QPen()
        pen.setColor(penColor)
        pen.setWidthF(self.lineWidth)
        painter.setPen(pen)

        try:
            startPoint=QtCore.QPointF(self.vertices[0][0]*300+self.width()/2.,self.vertices[0][2]*self.height()/2.)
        except IndexError:
            error_dialog = QtWidgets.QErrorMessage(self)
            error_dialog.showMessage('Please choose a data file first')
        else:
            def draw_whole_curve():
                start = 0
                end = self.numPt
                if self.isClosed == False:
                    start = 1
                for i in range(start, end):
                    p1 = QtCore.QPointF(self.vertices[i - 1][0] * 300 + self.width() / 2.,
                                        self.vertices[i - 1][2] * 300 + self.height() / 2.)
                    p2 = QtCore.QPointF(self.vertices[i][0] * 300 + self.width() / 2.,
                                        self.vertices[i][2] * 300 + self.height() / 2.)
                    painter.drawLine(p1, p2)

            def draw_composite_curve():
                if self.cut_points:
                    if self.isClosed == True:
                        for i in range(len(self.cut_points)):
                            penx = QtGui.QPen()
                            penx.setColor(self.color_list[i])
                            painter.setPen(penx)
                            row = self.vertices_matrix[i]
                            for j in range(len(row)-1):
                                p1 = QtCore.QPointF(row[j][0] * 300 + self.width() / 2.,
                                                    row[j][2] * 300 + self.height() / 2.)
                                p2 = QtCore.QPointF(row[j + 1][0] * 300 + self.width() / 2.,
                                                    row[j + 1][2] * 300 + self.height() / 2.)
                                painter.drawLine(p1, p2)

                    else:
                        for i in range(len(self.cut_points) + 1):
                            penx = QtGui.QPen()
                            penx.setColor(self.color_list[i])
                            painter.setPen(penx)
                            a = 0
                            b = self.numPt - 1
                            if i == 0:
                                b = self.cut_points[i]
                            elif i == len(self.cut_points):
                                a = self.cut_points[i - 1]
                            else:
                                a = self.cut_points[i - 1]
                                b = self.cut_points[i]
                            for j in range(a, b):
                                p1 = QtCore.QPointF(self.vertices[j][0] * 300 + self.width() / 2.,
                                                    self.vertices[j][2] * 300 + self.height() / 2.)
                                p2 = QtCore.QPointF(self.vertices[j + 1][0] * 300 + self.width() / 2.,
                                                    self.vertices[j + 1][2] * 300 + self.height() / 2.)
                                painter.drawLine(p1, p2)

                    painter.setPen(pen)
                else:
                    draw_whole_curve()


            if self.cut_points:
                draw_composite_curve()

            elif self.single_piece_mode==True:
                draw_whole_curve()

            elif self.fragment_mode==True:
                if self.start_index<self.end_index:
                    for i in range(self.start_index+1,self.end_index+1):   
                        p1=QtCore.QPointF(self.vertices[i-1][0]*300+self.width()/2.,self.vertices[i-1][2]*300+self.height()/2.)
                        p2=QtCore.QPointF(self.vertices[i][0]*300+self.width()/2.,self.vertices[i][2]*300+self.height()/2.)
                        painter.drawLine(p1,p2) 
                else:        
                    for i in range(self.start_index+1, self.numPt):
                        p1=QtCore.QPointF(self.vertices[i-1][0]*300+self.width()/2.,self.vertices[i-1][2]*300+self.height()/2.)
                        p2=QtCore.QPointF(self.vertices[i][0]*300+self.width()/2.,self.vertices[i][2]*300+self.height()/2.)
                        painter.drawLine(p1,p2)   
                    for i in range(self.end_index+1):   
                        p1=QtCore.QPointF(self.vertices[i-1][0]*300+self.width()/2.,self.vertices[i-1][2]*300+self.height()/2.)
                        p2=QtCore.QPointF(self.vertices[i][0]*300+self.width()/2.,self.vertices[i][2]*300+self.height()/2.)
                        painter.drawLine(p1,p2)

            if self.activateTangent == True:
                t = QtGui.QVector2D(self.tangents[i][0], self.tangents[i][1]).normalized() * 20.0
                p3 = QtCore.QPointF(t.x() + p2.x(), t.y() + p2.y())
                pen1 = QtGui.QPen()
                pen1.setColor(QtCore.Qt.blue)
                painter.setPen(pen1)
                painter.drawLine(p2, p3)
            if self.activateCurvature == True:
                normal = self.normals[i]
                p4 = QtCore.QPointF(normal.x() * self.curvatures[i] + p2.x(),
                                    normal.y() * self.curvatures[i] + p2.y())
                pen2 = QtGui.QPen()
                pen2.setColor(QtCore.Qt.green)
                painter.setPen(pen2)
                painter.drawLine(p2, p4)
            if self.activateCutPoint == True:  # just candidate cut points
                od = curve_fitting.sortCurvature(self.curvatures)
                # show first 5 cut point candidates based on curvatures on original curve
                for i in range(20):
                    index = od[-i][0]
                    p = self.vertices[index]
                    painter.drawText(p[0] * 300 + self.width() / 2., p[2] * 300 + self.height() / 2., str(index))
            painter.setPen(pen)

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
                angle=self.angles[i]
                a0u=a0u+(self.vertices[i][0]-self.center[0])*math.sin(angle)
                a0b=a0b+math.sin(angle)*math.sin(angle)
                b0u=b0u+(self.vertices[i][2]-self.center[2])*math.cos(angle)
                b0b=b0b+math.cos(angle)*math.cos(angle)
                
        except IndexError:
            error_dialog = QtWidgets.QErrorMessage(self)
            error_dialog.showMessage('Please choose a data file first')
        else:
            width=a0u/a0b
            height=b0u/b0b
            painter.drawEllipse(QtCore.QPointF(self.width()/2., self.height()/2.), width*300, height*300)


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

