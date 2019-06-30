# cross_section_transfer_script.py
#
# DESCRIPTION:
#    
#   This scripts add an attribute "transferMatrix" to every joint of target mesh, indicating
#   the transform matrix a cross-section curve should multiply with such that the curve
#   could align with the new mesh's surface.
#
#   For example, let's say we want the cross-section curves of 'source_male' to tansfer to the 'target_male',
#   we first select the group of cross-section curves that you want to transfer, i.e. 'source_cross_section_group'.
#   And also select the mesh you want the curves to transfer to, i.e. 'target_male'.
#   Then run this script, there will be a new group of cross-section curves created, named 'transfered_cross_section_group'.
#   In the viewport, the 'transfered_cross_section_group' will placed on 'target_male_mesh'.
# 
#   The output is the new group of cross section curves on target mesh
#   
#   Warning: this script only works for rest pose. Animation on skeleton will be ignored.
#------------------------------------------------- 

import maya.api.OpenMaya as om
import maya.cmds as cmds
import math
import sys

def maya_useNewAPI():
    pass
    
def getRotationMatrix(v0,v1):
    v0=v0.normalize()
    v1=v1.normalize()
    rotationAxis=v0^v1 # to determine a ratation axis
    theta=math.acos(v0*v1) # to find rotation angle
    # to build a quaternion
    s=math.sin(theta/2.0)
    qx=rotationAxis.x*s
    qy=rotationAxis.y*s
    qz=rotationAxis.z*s
    qw=math.cos(theta/2.0)
    # the transformation matrix is the quaternion as a 3 by 3
    # note that Matrix in maya is column-major
    rotMatrix=om.MMatrix()   
    rotMatrix.setElement(0,0,1-2*qy*qy-2*qz*qz)
    rotMatrix.setElement(1,0,2*qx*qy-2*qz*qw)
    rotMatrix.setElement(2,0,2*qx*qz+2*qy*qw)
    rotMatrix.setElement(3,0,0)
    rotMatrix.setElement(0,1,2*qx*qy+2*qz*qw)
    rotMatrix.setElement(1,1,1-2*qx*qx-2*qz*qz)
    rotMatrix.setElement(2,1,2*qy*qz-2*qx*qw)
    rotMatrix.setElement(3,1,0)
    rotMatrix.setElement(0,2,2*qx*qz-2*qy*qw)
    rotMatrix.setElement(1,2,2*qy*qz+2*qx*qw)
    rotMatrix.setElement(2,2,1-2*qx*qx-2*qy*qy)
    rotMatrix.setElement(3,2,0)
    rotMatrix.setElement(0,3,0)
    rotMatrix.setElement(1,3,0)
    rotMatrix.setElement(2,3,0)
    rotMatrix.setElement(3,3,1)
    
    return rotMatrix
       

selection=om.MSelectionList()
selection.add("Target_Belly")#belly joint is the root joint
rootJointObj=om.MObject()
rootJointDag=om.MDagPath()

MayaVersion=cmds.about(version=True)
if MayaVersion=="2017":
    selection.getDependNode(0,rootJointObj)
    selection.getDagPath(0,rootJointDag)
elif MayaVersion=="2018":
    rootJointObj=selection.getDependNode(0)
    rootJointDag=selection.getDagPath(0)
    
dagIter=om.MItDag()
dagIter.reset(rootJointDag,om.MItDag.kDepthFirst,om.MFn.kJoint)
transferMat=om.MMatrix()


while not dagIter.isDone():
    # Obtain the current item.
    curJointObj=dagIter.currentItem()
    # Make our MFnDagNode function set operate on the current DAG object.
    dagFn=om.MFnDagNode(curJointObj)
    #dgNodeFn=om.MFnDependencyNode(curJointObj)
    curJointName=dagFn.name()
    fatherCount=dagFn.parentCount()
    transferMat=om.MMatrix()
    
    # Add transfer matrix attribute to the joint, if it doesn't exist
    if not dagFn.hasAttribute("transferMatrix"):
        mttrFn=om.MFnMatrixAttribute()
        mAttr=mttrFn.create("transferMatrix","tfm",om.MFnMatrixAttribute.kDouble)
        mttrFn.readable=True
        mttrFn.storable=True # fairly consistent, won't change in compute()
        mttrFn.writable=False 
        dagFn.addAttribute(mAttr)
    
    # The root joint, a.k.a belly joint
    if fatherCount==0:
        # Get joint position
        # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
        jointPosition=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        # Transfer Matrix should indicate the translation
        translateMat=om.MMatrix.identity()
        translateMat.setElement(0,3,jointPosition[0])
        translateMat.setElement(1,3,jointPosition[1])
        translateMat.setElement(2,3,jointPosition[2])
        
        mPlug=dgNodeFn.findPlug('transferMatrix',False)
        sourceValueAsMObject = om.MFnMatrixData().create(translateMat)
        mPlug.setMObject( sourceValueAsMObject )
    # Require the parent's transfer Matrix attribute    
    if fatherCount==1:
        fatherJointObj=dagNodeFn.parent(0)
        fatherFn=om.MFnDagNode(fatherJointObj)
        fatherJointName=fatherFn.name()
        # Find father's transfer matrix attribute
        #fatherMatrixPlug=fatherFn.findPlug('transferMatrix',False)
        fatherTransferMatrix=fatherFn.attribute("transferMatrix")
        print fatherTransderMatrix
        p1=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        p0=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)        
        mPlug=dagNodeFn.findPlug('transferMatrix',False)
        translateVec=om.MVector((p1[0]-p0[0],p1[1]-p0[1],p1[2]-p1[2]))
        
        sourceValueAsMObject = om.MFnMatrixData().create(translateMat)
        mPlug.setMObject( sourceValueAsMObject )
        
    elif fatherCount>1:
        raise ValueError(curJointName+"has more than one parent")        
        
    dagIter.next()


nextBoneObj=rootBoneFn.child(0)
nextBoneFn=om.MFnDagNode(nextBoneObj)
nextBonePosition=cmds.xform(nextBoneFn.name(),absolute=True,query=True,worldSpace=True,rotatePivot=True)
nextBonePosition=om.MVector(nextBonePosition[0],nextBonePosition[1],nextBonePosition[2])
trgBoneVec=nextBonePosition-bonePosition

rm=getRotationMatrix(om.MVector((0,1,0)),trgBoneVec)

selection.clear()
selection.add('Source_Belly_cross_section_u_at_18_percentage1')
duplicateObj=om.MObject()
duplicateDag=om.MDagPath()
duplicateDag=selection.getDagPath(0)
transFn=om.MFnTransform(duplicateDag)
translateVector=om.MVector((2,2,2))
transFn.setTranslation(translateVector,om.MSpace.kWorld)