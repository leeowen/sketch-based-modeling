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
        qPlug=jointFn.findPlug('quaternionVec',False)
        """
        tPlug=jointFn.findPlug('translateVec',False)
        translateVec=om.MVector()
        translateVec.x=tPlug.child(0).asDouble()
        translateVec.y=tPlug.child(1).asDouble()            
        translateVec.z=tPlug.child(2).asDouble()
        """
        quaternionVec=om.MQuaternion()
        quarternionVec.x=qPlug.child(0).asDouble()
        quarternionVec.y=qPlug.child(1).asDouble()
        quarternionVec.z=qPlug.child(2).asDouble()
        quarternionVec.w=qPlug.child(3).asDouble()                        
        
        transformFn.rotateBy(quaternionVec,om.MSpace.kObject)
        #transformFn.translateBy(translateVec,om.MSpace.kWorld)
        return 
    else:
        if jointFn.hasParent():
            fatherJointObj=jointFn.parent(0)
            recursion(fatherJointObj,transformFn)
            
            qPlug=jointFn.findPlug('quaternionVec',False)
            quaternionVec=om.MQuaternion()
            quarternionVec.x=qPlug.child(0).asDouble()
            quarternionVec.y=qPlug.child(1).asDouble()
            quarternionVec.z=qPlug.child(2).asDouble()
            quarternionVec.w=qPlug.child(3).asDouble()   
            
            transformFn.rotateBy(quaternionVec,om.MSpace.kObject)
        else:
            pass
            

def linear_interpolate_3D(p1,p2,t):
    p1=om.MVector(p1[0],p1[1],p1[2])
    p2=om.MVector(p2[0],p2[1],p2[2])
    p=t*p2+(1.0-t)*p1
    return p


def getPositionAtU(jointObj,u):
    pos=om.MVector()
    jointFn=om.MFnDagNode(jointObj)
    # GET BONE DATA
    flag=1
    jointFn=om.MFnDagNode(jointObj)
    if jointFn.childCount()>0:
        next_bone_obj=jointFn.child(0)
        flag=-1
    elif jointFn.parentCount():
        next_bone_obj=jointFn.parent(0)
        flag=1
        
    nextjointFn=om.MFnDagNode(next_bone_obj)
    self.nextBone_name=nextjointFn.fullPathName().split('|')
    self.nextBone_name=self.nextBone_name[-1]
    self.nextBone_dagPath = nextBoneFn.getPath()
    
    # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
    bone_position=cmds.xform(self.bone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    bone_position=om.MVector(bone_position[0],bone_position[1],bone_position[2])
    nextBone_position=cmds.xform(self.nextBone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
    nextBone_position=om.MVector(nextBone_position[0],nextBone_position[1],nextBone_position[2])
    pos=linear_interpolate_3D(bone_position,nextBone_position,self.u_parameter)
    return pos
            

selection=om.MSelectionList()
selection.add("Target_Belly")#belly joint is the root joint
rootJointObj=om.MObject()
rootJointDag=om.MDagPath()

rootJointObj=selection.getDependNode(0)
rootJointDag=selection.getDagPath(0)
    
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
   
    # Add translate vector attribute to the joint, if it doesn't exist
    if not dagFn.hasAttribute("translateVec"):
        """
        # Make the translation vec an atrribute 
        nAttrFn = om.MFnUnitAttribute()
        cAttrFn = om.MFnCompoundAttribute()
        translateAttr = cAttrFn.create('translateVec', 'tsl')
        translateAttrX =nAttrFn.create('translateVecX', 'tlX', om.MFnNumericData.kDistance)
        translateAttrY =nAttrFn.create('translateVecY', 'tlY', om.MFnNumericData.kDistance)
        translateAttrZ =nAttrFn.create('translateVecZ', 'tlZ', om.MFnNumericData.kDistance)
        cAttrFn.setArray(True)
        cAttrFn.addChild(translateAttrX)
        cAttrFn.addChild(translateAttrY)
        cAttrFn.addChild(translateAttrZ)
        cAttrFn.setKeyable(True)
        dgFn.addAttribute(translateAttr)
        """
        
    if not dagFn.hasAttribute("quaternionVec"):
        uAttrFn = om.MFnNumericAttribute()
        qAttrFn = om.MFnCoumpoundAttribute()
        qAttr=qAttrFn.create('quaternionVec','qn')
        qAttrX=uAttrFn.create('quaternionVecX','qnx', om.MFnNumericData.kDouble)
        qAttrY=uAttrFn.create('quaternionVecY','qny', om.MFnNumericData.kDouble)
        qAttrZ=uAttrFn.create('quaternionVecZ','qnz', om.MFnNumericData.kDouble)
        qAttrW=uAttrFn.create('quaternionVecW','qnw', om.MFnNumericData.kDouble)
        qAttrFn.setArray(True)
        qAttrFn.addChild(qAttrX)
        qAttrFn.addChild(qAttrY)
        qAttrFn.addChild(qAttrZ)
        qAttrFn.addChild(qAttrW) 
        qAttrFn.writable=True  
        dgFn.addAttribute(qAttr)     
           
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
                quaternion=om.MVector.rotateTo(yAxis,newUp)              
    
    # If it is not root joint, a.k.a not Belly joint
    else:
        # Require the parent's transfer Matrix attribute 
        fatherJointObj=dagFn.parent(0)
        fatherFn=om.MFnDagNode(fatherJointObj)
        fatherJointName=fatherFn.name()

        # Find father joint's location to construct quaternion matrix and negative translate matrix      
        p1=cmds.xform(curJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        p0=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)        
        """
        mPlug=dagFn.findPlug('translateVec',False)
        translateVec.x=p0[0]
        translateVec.y=p0[1]
        translateVec.z=p0[2]
        """
        oldUp=om.MVector((p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]))
        
        
        # Find child joint's location to construct quaternion matrix and translate matrix
        if childCount==1:
            childJointObj=dagFn.child(0)
            childFn=om.MFnDagNode(childJointObj)
            childJointName=childFn.name()
            p2=cmds.xform(childJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
            newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
            quaternion=oldUp.rotateTo(newUp)
            
        # If has no child joint,i.e. an end joint, then no rotaion or translate Matrix is needed,
        # However, Hip joint is an exception
        elif childCount==0:
            if "Hip" in curJointName:
                p2=cmds.xform(fatherJointName,absolute=True,query=True,worldSpace=True,rotatePivot=True)
                newUp=om.MVector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
                quaternion=oldUp.rotateTo(newUp)
            else:
                pass
                
        # Have multiple children
        else:
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
    """
    cPlug=dgFn.findPlug('translateVec',False)
    for i in range(3):
        childPlug=cPlug.child(i)
        value=om.MDistance(translateVec[i],om.MDistance.kInvalid)
        childPlug.setMDistance(value)
    """
          
    mPlug=dagFn.findPlug('quaternionVec',False)        
    for i in range(4):
        childPlug=mPlug.child(i)
        value=quaternion[i]
        childPlug.setDouble(value)
    
    dagIter.next()


selection.clear()
originName='Source_Belly_cross_section_u_at_70_percentage'
selection.add(originName)
originDag=selection.getDagPath(0)
originFn=om.MFnDagNode(originDag)
"""
origPlug=om.MPlug()
origPlug=originFn.findPlug("objToWorld",False)
origMatObj=origPlug.asMObject()
origMat=om.MFnMatrixData(origMatObj).matrix()
worldToLocalMat=origMat.inverse()
dagModifier=om.MDagModifier()
"""
metaName=originName+'_meta'
selection.add(metaName)
metaDag=selection.getDagPath(1)
metaFn=om.MFnDagNode(metaDag)

metaTransformFn=om.MFnTransform(metaDag)
# The transform order is rotate(object space)-->translate

# Find the corresponding joint
jointName=metaName.split('_cross_section_')
selection.add(jointName)
jointObj=selection.getDependNode(2)
transformFn=om.MFnTransform(metaDag)
u=metaName.split('_')
u=int(u[6])
recursion(jointObj,transformFn,u)
positionAtU=getPositionAtU(jointObj,u)
transformFn.translateTo(positionAtU)

            