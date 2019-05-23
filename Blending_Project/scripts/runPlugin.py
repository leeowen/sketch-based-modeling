import maya.cmds as cmds
import sys

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")
cmds.crossSectionExtract(mu=0.5)
cmds.flushUndo()
cmds.scriptEditorInfo(clearHistory=True)