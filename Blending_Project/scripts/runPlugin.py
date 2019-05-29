import maya.cmds as cmds
import sys

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy, sympy


dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")
cmds.crossSectionExtract(mu=0.5,mmn="Source_male_meshShape",mbn="Source_LeftUpLeg")
#cmds.flushUndo()
cmds.crossSectionExtract(mu=0.5)
#cmds.scriptEditorInfo(clearHistory=True)