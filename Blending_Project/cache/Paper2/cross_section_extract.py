import maya.cmds as cmds
import maya.api.OpenMaya as OpenMaya
    
    
def VertexPositionOfEdge():
    #___________Selection___________
    # query the currently selected object
    selectionList=OpenMaya.MGlobal.getActiveSelectionList()	
    if selectionList.isEmpty():
        cmds.error("select edge from mesh please")
    iterObjects = OpenMaya.MItSelectionList(selectionList)

    edgeVertexSet=set()             
    while not iterObjects.isDone():   
        if not iterObjects.hasComponents():
            cmds.error("select edge from mesh please")
  
        dagPath,component=iterObjects.getComponent() 
        meshIter=OpenMaya.MItMeshEdge(dagPath,component)
        edgeVertexList=[]
        
        #___________Query each vertex position ___________
        while not meshIter.isDone():
            pt=meshIter.point(0,OpenMaya.MSpace.kWorld)# use the world space
            pt=(pt.x,pt.y,pt.z)
            edgeVertexList.append(pt)
            pt=meshIter.point(1,OpenMaya.MSpace.kWorld)
            pt=(pt.x,pt.y,pt.z)
            edgeVertexList.append(pt)
            meshIter.next()
        
        edgeVertexSet=set(edgeVertexList)
        iterObjects.next()

    return edgeVertexSet
    
 
VertexPositionOfEdge()