import maya.api.OpenMaya as om

selection=om.MGlobal.getActiveSelectionList()
meshIter=om.MItSelectionList(selection)

dagPath,component=meshIter.getComponent()
print dagPath,component

