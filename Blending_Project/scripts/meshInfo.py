import maya.api.OpenMaya as om

selection=om.MGlobal.getActiveSelectionList()
meshIter=om.MItSelectionList(selection)
dagPath,component=meshIter.getComponent()
vertIter=om.MItMeshVertex(dagPath,component)

while not vertIter.isDone():
    print vertIter.position(om.MSpace.kObject)
    vertIter.next()