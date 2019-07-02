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
#   Warning: this script only works for rest pose. Animations on skeleton will be ignored.
#------------------------------------------------- 

import maya.api.OpenMaya as om
import maya.cmds as cmds
import math
import sys

def maya_useNewAPI():
    pass
    
def getTranslateMatrix(p):
    tMat=om.MMatrix()
    tMat.setElement(0,0,1)
    tMat.setElement(0,1,0)
    tMat.setElement(0,2,0)
    tMat.setElement(0,3,p[0])
    tMat.setElement(1,0,0)
    tMat.setElement(1,1,1)
    tMat.setElement(1,2,0)
    tMat.setElement(1,3,p[1])
    tMat.setElement(2,0,0)
    tMat.setElement(2,1,0)
    tMat.setElement(2,2,1)
    tMat.setElement(2,3,p[2])
    tMat.setElement(3,0,0)
    tMat.setElement(3,1,0)
    tMat.setElement(3,2,0)
    tMat.setElement(3,3,1)
    return tMat
    
def getQuaternionMatrix(v0,v1):
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
       
def decomposeTransformationMatrix(mat):
    translateVec=om.MVector((mat.getElement(3,0),mat.getElement(3,1),mat.getElement(3,2)))

    sx=om.MVector((mat.getElement(0,0),mat.getElement(0,1),mat.getElement(0,2))).length()
    sx=om.MVector((mat.getElement(1,0),mat.getElement(1,1),mat.getElement(1,2))).length()
    sx=om.MVector((mat.getElement(2,0),mat.getElement(2,1),mat.getElement(2,2))).length()
    scaleVec=om.MVector(sx,sy,sz)

    rotationMat=om.MMatrix()
    rotationMat.setElement(0,0,mat.getElement(0,0)/sx)
    rotationMat.setElement(1,0,mat.getElement(1,0)/sx)
    rotationMat.setElement(2,0,mat.getElement(2,0)/sx)
    rotationMat.setElement(3,0,0)
    rotationMat.setElement(0,1,mat.getElement(0,1)/sy)
    rotationMat.setElement(1,1,mat.getElement(1,1)/sy)
    rotationMat.setElement(2,1,mat.getElement(2,1)/sy)
    rotationMat.setElement(3,1,0)
    rotationMat.setElement(0,2,mat.getElement(0,2)/sz)
    rotationMat.setElement(1,2,mat.getElement(1,2)/sz)
    rotationMat.setElement(2,2,mat.getElement(2,2)/sz)
    rotationMat.setElement(3,2,0)
    rotationMat.setElement(0,3,0)
    rotationMat.setElement(1,3,0)
    rotationMat.setElement(2,3,0)
    rotationMat.setElement(3,3,1)
    
    return translateVec,scaleVec,rotationMat
    

selection=om.MSelectionList()
selection.add("Target_Belly")#belly joint is the root joint
rootJointObj=om.MObject()
rootJointDag=om.MDagPath()

"""
MayaVersion=cmds.about(version=True)
if MayaVersion=="2017":
    selection.getDependNode(0,rootJointObj)
    selection.getDagPath(0,rootJointDag)
elif MayaVersion=="2018":
    rootJointObj=selection.getDependNode(0)
    rootJointDag=selection.getDagPath(0)
"""
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
    childCount=dagFn.childCount()
    transferMat=om.MMatrix().setToIdentity()
    translateMat=om.MMatrix().setToIdentity()
    quarternionMat=om.MMatrix().setToIdentity()
    translateBackMat=om.MMatrix().setToIdentity()
    
    # Add transfer matrix attribute to the joint, if it doesn't exist
    if not dagFn.hasAttribute("transferMatrix"):
        mttrFn=om.MFnMatrixAttribute()
        mAttr=mttrFn.create("transferMatrix","tfm",om.MFnMatrixAttribute.kDouble)
        mttrFn.readable=True
        mttrFn.storable=True # fairly consistent, won't change in compute()
        mttrFn.writable=False 
        dagFn.addAttribute(mAttr)
    
    # The root joint, a.k.a belly joint
    if "Belly" in curJointName:
        # Get joint position
        # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
        p0=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        
        for i in range(childCount):
            childJointObj=dagFn.child(i)
            childFn=om.MFnDagNode(childJointObj)
            childJointName=childFn.name()
            if "Chest" in childJointName: 
                p1=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                quaternionMat=getQuaternionMatrix(om.MVector((0,1,0)),newUp)
            
        # Transfer Matrix should indicate the translation
        translateMat.setElement(0,3,p0[0])
        translateMat.setElement(1,3,p0[1])
        translateMat.setElement(2,3,p0[2])
        
        # The transform order is rotate-->translate
        # The transform matrix M is after the point/vector P (P' = P x M)" 
        transferMat=quarternionMat*translateMat   
        
        # Networked plugs contain the actual connection data between two nodes, 
        # and maybe some other flags as required, such as whether the plug is locked. 
        # Maya is free to create and destroy them as it needs to, and no notifications are sent when that happens. 
        # Non-networked plugs can be thought of as named references to a plug. They belong to you and have the lifespan you dictate.
        # In C++ terms a non-networked plug is like a copy of an object and a networked plug is like a pointer to it. 
        mPlug=dagFn.findPlug('transferMatrix',False)
        sourceValueAsMObject = om.MFnMatrixData().create(transferMat)
        mPlug.setMObject( sourceValueAsMObject )   
    
    # If it is not root joint, a.k.a not Belly joint
    else:
        # Require the parent's transfer Matrix attribute 
        fatherJointObj=dagFn.parent(0)
        fatherFn=om.MFnDagNode(fatherJointObj)
        fatherJointName=fatherFn.name()

        fatherPlug=fatherFn.findPlug('transferMatrix',False)
        fatherMatObj=fatherPlug.asMObject()
        fatherMat=om.MFnMatrixData(fatherMatObj).matrix()

        # Find father joint's location to construct quaternion matrix and negative translate matrix      
        p1=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        p0=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)        
        mPlug=dagFn.findPlug('transferMatrix',False)
        translateBackVec=om.MVector((-p0[0],-p0[1],-p0[2]))
        translateBackMatrix=getTranslateMatrix(translateBackVec)
        oldUp=om.MVector((p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]))
        
        # Find child joint's location to construct quaternion matrix and translate matrix
        if childCount==1:
            childJointObj=dagFn.child(0)
            childFn=om.MFnDagNode(childJointObj)
            childJointName=childFn.name()
            p2=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
            newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
            quaternionMat=getQuaternionMatrix(oldUp,newUp)
            
        # If has no child joint,i.e. an end joint, then no rotaion or translate Matrix is needed,
        # However, Hip joint is an exception
        elif childCount==0:
            if "Hip" in curJointName:
                p2=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
                quaternionMat=getQuaternionMatrix(oldUp,newUp)
            else:
                quarternionMat.setToIdentity()
                
        # Have multiple children
        else:
            for i in range(childCount):
                childJointObj=dagFn.child(i)
                childFn=om.MFnDagNode(childJointObj)
                childJointName=childFn.name()
                if "Head" in childJointName: 
                    p2=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                    newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                    quaternionMat=getQuaternionMatrix(om.MVector((0,1,0)),newUp)                                        
        
        translateMat=getTranslateMatrix(p1)        
        transferMat=fatherMat*translateBackMat*quaternionMat*translateMat    
        sourceValueAsMObject = om.MFnMatrixData().create(transferMat)
        mPlug.setMObject( sourceValueAsMObject )        
        
    dagIter.next()


selection.clear()
originName='Source_Belly_cross_section_u_at_70_percentage'
selection.add(originName)
originDag=selection.getDagPath(0)
originFn=om.MFnDagNode(originDag)
origPlug=om.MPlug()
origPlug=originFn.findPlug("objToWorld",False)
origMatObj=origPlug.asMObject()
origMat=om.MFnMatrixData(origMatObj).matrix()
worldToLocalMat=origMat.inverse()
dagModifier=om.MDagModifier()
instanceObj=dagModifier.createNode( 'transform' )
name=originName.split("Source")
duplicateName="ideal"+name[1]
dagModifier.renameNode(instanceObj,duplicateName)
print worldToLocalMat.getElement(0,1)
translationVec,scaleVec,rotationMat=decomposeTransformationMatrix(worldToLocalMat)
#duplicateMatrix=om.MTransformationMatrix(worldToLocalMat)
#translation=duplicateMatrix.translation(om.MSpace.kWorld)

#div=40
#eps=om.MPointsArray()
#for i in range (div):
#    u=1.0/div*i
#    pt=origFn.getPointAtParam(u,pm.MSpace.kWorld)
#    eps.append()

