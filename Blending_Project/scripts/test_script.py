import maya.api.OpenMaya as om
import maya.cmds as cmds
from PySide2 import QtCore
from PySide2 import QtWidgets,
from PySide2 import wrapInstance
import maya.api.OpenMayaUI as omui
import math,sys
sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy


def maya_useNewAPI():
    pass


cmds.flushUndo()
cmds.scriptEditorInfo(clearHistory=True)