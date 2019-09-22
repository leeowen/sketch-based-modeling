import maya.cmds as cmds
from PySide2 import QtCore, QtGui, QtWidgets
from shiboken2 import wrapInstance
import maya.OpenMayaUI as omui

    
def maya_useNewAPI():
    pass


def maya_main_window():
    main_window_ptr=omui.MQtUtil.mainWindow()
    return wrapInstance(long(main_window_ptr),QtWidgets.QWidget)
        

class ImportFileUI(QtWidgets.QDialog):
    
    def __init__(self,parent=maya_main_window()):
        super(ImportFileUI,self).__init__(parent)
            
        self.setObjectName("importFileDialogue")
        self.setWindowTitle("import cross-section contours")
        
        self.create_widgets()   
        self.create_layout()
        self.create_connection()
        
        
    def create_widgets(self):
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
        self.confirm_button=QtWidgets.QPushButton('Confirm')
        self.cancel_button=QtWidgets.QPushButton('Cancel')
 

    def create_layout(self):
        main_layout=QtWidgets.QVBoxLayout(self)
        
        file_layout=QtWidgets.QHBoxLayout()
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.filePath_lineEdit)
        file_layout.addWidget(self.file_button)
        
        button_layout=QtWidgets.QHBoxLayout()
        button_layout.addWidget(self.confirm_button)
        button_layout.addWidget(self.cancel_button)
        
        main_layout.addLayout(file_layout)
        main_layout.addLayout(button_layout)
        
        
    def create_connection(self):
        self.file_button.clicked.connect(self.open_file_dialog)
        self.cancel_button.clicked.connect(self.close)
        self.confirm_button.clicked.connect(self.data_process)
        
    
    def open_file_dialog(self):
        dirPath=cmds.workspace(q=True, rootDirectory=True )
        FILE_FILTERS="data(*.dat);;All Files(*.*)"
        selected_filter="data(*.dat)"# default filter, also store last selected filter and can be used as the default filter for next select
        file_paths,selected_filter=QtWidgets.QFileDialog.getOpenFileNames(self, 'Select File',dirPath+'data',FILE_FILTERS,selected_filter)
        # check if user has cancel the dialog by checking if file_path is none
        if file_paths:
            #self.filePath_lineEdit.setText(file_paths)
            self.read_files(file_paths)
            
        self.close()
        

    def read_files(self,file_paths):
        self.numCr=0
        for file_path in file_paths:
            f=open(file_path,'r')
            content=f.readlines()
            
            for line in content:
                p=line.split()
               
            self.numCr++
            f.close()
        
        
if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("importFileDialogue",exists=True):
        cmds.deleteUI("importFileDialogue",wnd=True)
        
    w=ImportFileUI()
    w.show() 
