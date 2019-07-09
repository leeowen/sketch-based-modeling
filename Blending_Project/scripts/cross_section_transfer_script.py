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
   

def linear_interpolate_3D(p1,p2,t):
    p1=om.MVector(p1[0],p1[1],p1[2])
    p2=om.MVector(p2[0],p2[1],p2[2])
    p=t*p2+(1.0-t)*p1
    return p


def getPositionAtU(jointObj,u_parameter):
    pos=om.MVector()
    jointFn=om.MFnDagNode(jointObj)
    jointName=jointFn.name()
    nextJointObj=om.MObject()
    # GET BONE DATA
    jointFn=om.MFnDagNode(jointObj)
    if jointFn.childCount()>1:
        if 'Belly' in jointName:
            for i in range(jointFn.childCount()):
                childObj=jointFn.child(i)
                childName=om.MFnDependencyNode(childObj).name()
                if 'Chest' in childName:
                    nextJointObj=childObj
        elif 'Neck' in jointName:
            for i in range(jointFn.childCount()):
                childObj=jointFn.child(i)
                childName=om.MFnDependencyNode(childObj).name()
                if 'Head' in childName:
                    nextJointObj=childObj
    elif jointFn.childCount()==1:
        nextJointObj=jointFn.child(0)
    elif jointFn.childCount()==0 and jointFn.parentCount()==1:
        nextJointObj=jointFn.parent(0)
                
    nextJointFn=om.MFnDagNode(nextJointObj)
    nextBoneName=nextJointFn.name()
    nextBoneDagPath = nextJointFn.getPath()
    
    # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
    jointPosition=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    jointPosition=om.MVector(jointPosition[0],jointPosition[1],jointPosition[2])
    nextjointPosition=cmds.xform(nextBoneName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    nextjointPosition=om.MVector(nextjointPosition[0],nextjointPosition[1],nextjointPosition[2])
    pos=linear_interpolate_3D(jointPosition,nextjointPosition,u_parameter)
    positionAtU=om.MVector(pos[0],pos[1],pos[2])    
    return positionAtU
            

def getQuaternion(jointObj):
    jointFn=om.MFnDagNode(jointObj)
    jointName=jointFn.name()
    quaternion=om.MQuaternion()
    # The root joint, a.k.a belly joint
    if "Belly" in jointName:
        # Get joint position
        # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
        p0=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        
        for i in range(jointFn.childCount()):
            childJointObj=jointFn.child(i)
            childFn=om.MFnDagNode(childJointObj)
            childJointName=childFn.name()
            if "Chest" in childJointName: 
                p1=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                quaternion=om.MVector.kYaxisVector.rotateTo(newUp)   
                break           
    
    elif 'Hip' in jointName: 
        p1=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        fatherObj=jointFn.parent(0)
        fatherFn=om.MFnDependencyNode(fatherObj)
        fatherJointName=fatherFn.name()
        p0=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
        quaternion=om.MVector.kYnegAxisVector.rotateTo(newUp) 
    elif 'Neck' in jointName:
        p0=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        for i in range(jointFn.childCount()):
            childObj=jointFn.child(0)
            childJointName=om.MFnDependencyNode(childObj).name()
            if 'Head' in childJointName:
                p1=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                quaternion=om.MVector.kYaxisVector.rotateTo(newUp) 
                break
    elif jointFn.parentCount==1 and jointFn.childCount==1:
        p0=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        childObj=jointFn.parent(0)
        childFn=om.MFnDependencyNode(childObj)
        childJointName=childFn.name()
        p1=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
        if newUp.y>0:
            quaternion=om.MVector.kYaxisVector.rotateTo(newUp) 
        else:
            quaternion=om.MVector.kYnegAxisVector.rotateTo(newUp)
    return quaternion
    

selection.clear()
selection=om.MGlobal.getSelectionListByName("*_percentage_meta")

for i in range(selection.length()):
    metaCrv=selection.getDependNode(i)
    metaFn=om.MFnDagNode(metaCrv)
    metaName=metaFn.name()   
    # Find the corresponding joint
    strList=metaName.split('_')
    jointName="Target_"+strList[1]
    sList=om.MSelectionList()
    sList.add(jointName)
    jointObj=sList.getDependNode(0)
    u=int(strList[6])
    u=u/100.0
    qn=getQuaternion(jointObj)
    positionAtU=getPositionAtU(jointObj,u)
    transformFn=om.MFnTransform(metaCrv) 
    transformFn.setTranslation(positionAtU,om.MSpace.kObject)
    transformFn.rotateBy(qn,om.MSpace.kObject)





            