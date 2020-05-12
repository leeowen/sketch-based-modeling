import maya.cmds as cmds
cmds.flushUndo()
cmds.scriptEditorInfo(clearHistory=True)
dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.unloadPlugin(dirPath+"cross_section_extract_cmd.py")

reload(curve_fitting)
