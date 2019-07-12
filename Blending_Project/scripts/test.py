import maya.cmds as cmds
from PySide2 import QtCore, QtGui, QtWidgets
from shiboken2 import wrapInstance
import maya.OpenMayaUI as omui
import math,sys


class Canvas(QtWidgets.QWidget):
    backgroundColor=QtGui.QColor(100,122,200)
    def __init__(self,parent):
        super(Canvas,self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.WindowType.Window)
        self.resize(800,600)
        self.setObjectName('test')
        self._pixmap = QtGui.QPixmap("/home/i7664633/Ismail_color_scheme.jpg")
        
        self.show()
        
    def paintevent(self,evt):
        painter=QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing,True)
        # Draw the background
        #painter.setBrush(Canvas.backgroundColor)
        #painter.drawRect(self.rect())
        print self._pixmap.isNull()
        painter.drawPixmap(evt.rect(), self._pixmap)

    
def maya_main_window():
    main_window_ptr=omui.MQtUtil.mainWindow()
    return wrapInstance(long(main_window_ptr),QtWidgets.QWidget)
    
    
if __name__=="__main__":
    # Check to see if the UI already exists and if so, delete
    if cmds.window("test",exists=True):
       cmds.deleteUI("test",wnd=True)
    """
    try:
        w.close()
        w.deleteLater()
    except:
        pass
    """
        
        
    w=Canvas(maya_main_window())
    w.show()  
    w.raise_() 
