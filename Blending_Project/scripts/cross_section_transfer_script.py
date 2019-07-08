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
   
def recursion(jointObj,transformFn,u):
    jointFn=om.MFnDagNode(jointObj)
    jointName=jointFn.name()
    if 'Belly' in jointName or 'Hip' in jointName:
        if jointFn.hasAttribute('quaternionVec'):
            qPlug=jointFn.findPlug('quaternionVec',False)
            quaternionVec=om.MQuaternion()
            quaternionVec.x=qPlug.child(0).asDouble()
            quaternionVec.y=qPlug.child(1).asDouble()
            quaternionVec.z=qPlug.child(2).asDouble()
            quaternionVec.w=qPlug.child(3).asDouble() 
        else:
            raise ValueError( "can't find quaternion attribute of node "+jointName )                      
        
        transformFn.rotateBy(quaternionVec,om.MSpace.kObject)
        return 
    else:
        if jointFn.hasParent():
            fatherJointObj=jointFn.parent(0)
            recursion(fatherJointObj,transformFn)
            if jointFn.hasAttribute('quaternionVec'):
                qPlug=jointFn.findPlug('quaternionVec',False)
                quaternionVec=om.MQuaternion()
                quarternionVec.x=qPlug.child(0).asDouble()
                quarternionVec.y=qPlug.child(1).asDouble()
                quarternionVec.z=qPlug.child(2).asDouble()
                quarternionVec.w=qPlug.child(3).asDouble()   
            else:
                print "can't find quaternion attribute"
            
            transformFn.rotateBy(quaternionVec,om.MSpace.kObject)
        else:
            pass
            

def linear_interpolate_3D(p1,p2,t):
    p1=om.MVector(p1[0],p1[1],p1[2])
    p2=om.MVector(p2[0],p2[1],p2[2])
    p=t*p2+(1.0-t)*p1
    return p


def getPositionAtU(jointObj,u_parameter):
    pos=om.MVector()
    jointFn=om.MFnDagNode(jointObj)
    # GET BONE DATA
    flag=1
    jointFn=om.MFnDagNode(jointObj)
    if jointFn.childCount()>0:
        nextJointObj=jointFn.child(0)
        flag=-1
    elif jointFn.parentCount():
        nextJointObj=jointFn.parent(0)
        flag=1
        
    nextJointFn=om.MFnDagNode(nextJointObj)
    nextBoneName=nextJointFn.name()
    nextBoneDagPath = nextJointFn.getPath()
    
    # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
    jointPosition=cmds.xform(jointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    jointPosition=om.MVector(jointPosition[0],jointPosition[1],jointPosition[2])
    nextjointPosition=cmds.xform(nextBoneName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    nextjointPosition=om.MVector(nextjointPosition[0],nextjointPosition[1],nextjointPosition[2])
    pos=linear_interpolate_3D(jointPosition,nextjointPosition,u_parameter)
    return pos
            

def processJoint(rootJointDag):
    dagIter=om.MItDag()
    dagIter.reset(rootJointDag,om.MItDag.kDepthFirst,om.MFn.kJoint)

    while not dagIter.isDone():
        # Obtain the current item.
        curJointObj=dagIter.currentItem()
        # Make our MFnDagNode function set operate on the current DAG object.
        dagFn=om.MFnDagNode(curJointObj)
        #dgNodeFn=om.MFnDependencyNode(curJointObj)
        curJointName=dagFn.name()
        fatherCount=dagFn.parentCount()
        childCount=dagFn.childCount()
        quaternion=om.MQuaternion()
            
        if not dagFn.hasAttribute("quaternionVec"):
            uAttrFn = om.MFnNumericAttribute()
            qAttrFn = om.MFnCompoundAttribute()
            qAttr=qAttrFn.create('quaternionVec','qn')
            qAttrX=uAttrFn.create('quaternionVecX','qnx', om.MFnNumericData.kDouble)
            qAttrY=uAttrFn.create('quaternionVecY','qny', om.MFnNumericData.kDouble)
            qAttrZ=uAttrFn.create('quaternionVecZ','qnz', om.MFnNumericData.kDouble)
            qAttrW=uAttrFn.create('quaternionVecW','qnw', om.MFnNumericData.kDouble)
            qAttrFn.addChild(qAttrX)
            qAttrFn.addChild(qAttrY)
            qAttrFn.addChild(qAttrZ)
            qAttrFn.addChild(qAttrW) 
            qAttrFn.writable=True  
            qAttrFn.readable=True
            qAttrFn.keyable=True
            qAttrFn.storable=True
            qAttrFn.displayable=True
            qAttrFn.channelBox=True
            qAttrFn.dynamic=True
            dagFn.addAttribute(qAttr)   
     
        # The root joint, a.k.a belly joint
        if "Belly" in curJointName:
            # Get joint position
            # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
            p0=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
            #translateVec=om.MVector((p0[0],p0[1],p0[2]))
            
            for i in range(childCount):
                childJointObj=dagFn.child(i)
                childFn=om.MFnDagNode(childJointObj)
                childJointName=childFn.name()
                if "Chest" in childJointName: 
                    p1=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                    newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                    quaternion=om.MVector.kYaxisVector.rotateTo(newUp)              
        
        # If it is not root joint, a.k.a not Belly joint
        else:
            # Require the parent's transfer Matrix attribute 
            fatherJointObj=dagFn.parent(0)
            fatherFn=om.MFnDagNode(fatherJointObj)
            fatherJointName=fatherFn.name()
    
            # Find father joint's location to construct quaternion matrix and negative translate matrix      
            p1=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
            p0=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)        
            oldUp=om.MVector((p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]))
            
        
            # Find child joint's location to construct quaternion matrix and translate matrix
            if childCount==1:
                childJointObj=dagFn.child(0)
                childFn=om.MFnDagNode(childJointObj)
                childJointName=childFn.name()
                p2=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
                quaternion=oldUp.rotateTo(newUp)
                
            # If has no child joint,i.e. an end joint, then no rotaion is needed,
            # However, Hip joint is an exception
            elif childCount==0:               
                if "Hip" in curJointName:
                    p2=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                    newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
                    quaternion=om.MVector.kYaxisVector.rotateTo(newUp)                    
            # Have multiple children, e.g. Neck 
            else:
                if not "Neck" in curJointName:
                    raise ValueError(curJointName+"shouldn't have multiple children")
                for i in range(childCount):
                    childJointObj=dagFn.child(i)
                    childFn=om.MFnDagNode(childJointObj)
                    childJointName=childFn.name()
                    if "Head" in childJointName: 
                        p2=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                        newUp=om.MVector(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2])
                        quaternion=oldUp.rotateTo(newUp)                                            
        
        # Networked plugs contain the actual connection data between two nodes, 
        # and maybe some other flags as required, such as whether the plug is locked. 
        # Maya is free to create and destroy them as it needs to, and no notifications are sent when that happens. 
        # Non-networked plugs can be thought of as named references to a plug. They belong to you and have the lifespan you dictate.
        # In C++ terms a non-networked plug is like a copy of an object and a networked plug is like a pointer to it.                                     
        
        mPlug=dagFn.findPlug('quaternionVec',False)                
        for i in range(4):
            childPlug=mPlug.child(i)
            value=quaternion[i]
            childPlug.setDouble(value)
        
        dagIter.next()

selection=om.MSelectionList()
selection.add("Target_Belly")#belly joint is the root joint
rootJointObj=om.MObject()
rootJointDag=om.MDagPath()

rootJointObj=selection.getDependNode(0)
rootJointDag=selection.getDagPath(0)

processJoint(rootJointDag)    

selection.clear()
strList=cmds.select("source_meta_cross_section_group", hierarchy=True)
print strList
“”“
for originCrv in strList:
    originFn=om.MFnDagNode(originCrv)
    metaName=originFn.name()    
    # Find the corresponding joint
    strList=metaName.split('_')
    jointName="Target_"+strList[1]
    
    selection.add(jointName)
    jointObj=selection.getDependNode(2)
    transformFn=om.MFnTransform(metaDag)
    u=int(strList[6])
    u=u/100.0
    recursion(jointObj,transformFn,u)
    positionAtU=getPositionAtU(jointObj,u)
    
    transformFn.translateBy(positionAtU,om.MSpace.kWorld)
”“”
            