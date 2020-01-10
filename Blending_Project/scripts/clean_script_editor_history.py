import maya.cmds as cmds
#cmds.flushUndo()
cmds.scriptEditorInfo(clearHistory=True)
reload(curve_fitting)