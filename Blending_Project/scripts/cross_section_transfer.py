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
    """
    q=om.MQuaternion(theta,rotationAxis)
    return q
    """
    
class Tree(object):
    "Generic tree node."
    def __init__(self, name='root', children=None):
        self.name = name
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)
    def __repr__(self):
        return self.name
    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)
       

"""
skeleton=Tree('Belly',[Tree('Hip'),
                       Tree('LeftThigh',[Tree('LeftLeg',[Tree('LeftFoot',[Tree('LeftToe')])])],
                       Tree('RightThigh',[Tree('RightLeg',[Tree('RightFoot',[Tree('RightToe')])])],
                       Tree('Chest',[Tree('Neck'),[Tree('LeftArm'),
                                                   Tree('RightArm')])])))])        
"""

selection=om.MSelectionList()
selection.add("Target_Belly")#belly joint is the root joint
rootJointObj=om.MObject()
rootJointDag=om.MDagPath()
selection.getDependNode(0,rootJointObj)
selection.getDagPath(0,rootJointDag)

dagIter=om.MItDag()
dagIter.reset(rootJointDag,om.MItDag.kDepthFirst,om.MFn.kJoint)
dagNodeFn=om.MFnDagNode()

curJointObj=dagIter.currentItem()
dagNodeFn.setObject( curJointObj )
fatherCount=dagNodeFn.parentCount()
transferMat=om.MMatrix()
if fatherCount==0:
    jointPosition=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    # Transfer Matrix should indicate the translation
    translateMat=om.MMatrix.identity()
    translateMat.setElement(0,3,jointPosition[0])
    translateMat.setElement(1,3,jointPosition[1])
    translateMat.setElement(2,3,jointPosition[2])
    # Add transfer matrix attribute to the node
    mttrFn=om.MFnMatrixAttribute()
    mAttr=mttrFn.create("transferMatrix","tfm")
    mttrFn.readable=True
    mttrFn.storable=True
    mttrFn.writable=False
    dagNodeFn.addAttribute(mAttr)
    
    mPlug=dagNodeFn.findPlug('transferMatrix',False)
    sourceValueAsMObject = om.MFnMatrixData().create(translateMat)
    mPlug.setMObject( sourceValueAsMObject )
    
    

while not dagIter.isDone():
    # Obtain the current item.
    curJointObj=iter.currentItem()
    # Make our MFnDagNode function set operate on the current DAG object.
    dagNodeFn.setObject( curJointObj )
    curJointName=dagNodeFn.name()
    fatherCount=dagNodeFn.parentCount()
    transferMat=om.MMatrix()
    if fatherCount>0:
        # Require the parent's transfer Matrix attribute
        if fatherCount==1:
            fatherJointObj=dagNodeFn.parent(0)
            fatherFn=om.MDagNodeFn(fatherJointObj)
            # Find father's transfer matrix attribute
            fatherTransferMatrix=fatherFn.attribute("transferMatrix")
            
            # Add transfer matrix attribute to the node
            mttrFn=om.MFnMatrixAttribute()
            mAttr=mttrFn.create("transferMatrix","tfm")
            mttrFn.readable=True
            mttrFn.storable=True
            mttrFn.writable=False
            dagNodeFn.addAttribute(mAttr)
            
            mPlug=dagNodeFn.findPlug('transferMatrix',False)
            sourceValueAsMObject = om.MFnMatrixData().create(translateMat)
            mPlug.setMObject( sourceValueAsMObject )
        elif fatherCount==0:
            raise ValueError(curJointName+"has more than one parent")
    else:
        # It is the root joint, a.k.a belly joint
        # Get joint position
        # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
        jointPosition=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        # Transfer Matrix should indicate the translation
        translateMat=om.MMatrix.identity()
        translateMat.setElement(0,3,jointPosition[0])
        translateMat.setElement(1,3,jointPosition[1])
        translateMat.setElement(2,3,jointPosition[2])
        # Add transfer matrix attribute to the node
        mttrFn=om.MFnMatrixAttribute()
        mAttr=mttrFn.create("transferMatrix","tfm")
        mttrFn.readable=True
        mttrFn.storable=True
        mttrFn.writable=False
        dagNodeFn.addAttribute(mAttr)
        
        mPlug=dagNodeFn.findPlug('transferMatrix',False)
        sourceValueAsMObject = om.MFnMatrixData().create(translateMat)
        mPlug.setMObject( sourceValueAsMObject )
        
        
    dagIter.next()


rootBoneFn=om.MFnDagNode(rootJointDag)
bonePosition=cmds.xform(rootBoneFn.name(),absolute=True,query=True,worldSpace=True,rotatePivot=True)
bonePosition=om.MVector(bonePosition[0],bonePosition[1],bonePosition[2])
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