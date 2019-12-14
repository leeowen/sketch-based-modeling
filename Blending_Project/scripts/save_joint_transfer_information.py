import maya.api.OpenMaya as om
import pymel.core as pmc
import math
import sys


def maya_useNewAPI():
    pass


def get_child_joint(jointObj):
    jointFn = om.MFnDagNode(jointObj)
    jointName = jointFn.name()
    childJointObj = om.MObject()

    if jointFn.childCount() > 1:
        if 'Belly' in jointName:
            for i in range(jointFn.childCount()):
                childObj = jointFn.child(i)
                childName = om.MFnDependencyNode(childObj).name()
                if 'Chest' in childName:
                    childJointObj = childObj
        elif 'Neck' in jointName:
            for i in range(jointFn.childCount()):
                childObj = jointFn.child(i)
                childName = om.MFnDependencyNode(childObj).name()
                if 'Head' in childName:
                    childJointObj = childObj
    elif jointFn.childCount() == 1:
        childJointObj = jointFn.child(0)
    elif jointFn.childCount() == 0:
        return

    childJointFn = om.MFnDagNode(childJointObj)
    childJointName = childJointFn.name()
    childJointDagPath = childJointFn.getPath()
    return childJointName, childJointObj, childJointDagPath


def get_parent_joint(jointObj):
    jointFn = om.MFnDagNode(jointObj)
    parentJointObj = om.MObject()

    # The root joint has no parent, therefore, return right away
    if jointFn.parentCount == 0:
        return
    elif jointFn.parentCount == 1:
        parentJointObj = jointFn.parent(0)

    parentJointFn = om.MFnDependencyNode(fatherObj)
    parentJointName = parentJointFn.name()
    parentJointDagPath = parentJointFn.getPath()

    return parentJointName, parentJointObj, parentJointDagPath


def get_joint_transform(jointObj):
    jointFn = om.MFnDependencyNode(jointObj)
    jointName = jointFn.name()
    # get next joint
    if 'Hip' in jointName:
        nextJointName,nextJointObj,nextJointDagPath = get_parent_joint(jointObj)
    elif ('Tip' in jointName) or ('Toe' in jointName) or ('Hand' in jointName):
        raise ValueError('Foot toes, hands or tip of the head does not need to get transform information')
    else:
        nextJointName, nextJointObj, nextJointDagPath = get_child_joint(jointObj)

    p0 = cmds.xform(jointName, absolute=True, query=True, worldSpace=True, rotatePivot=True)
    p1 = cmds.xform(nextJointName, absolute=True, query=True, worldSpace=True, rotatePivot=True)
    quaternion = getQuaternion(p0,p1)
    return p0,p1,quaternion


def getQuaternion(p0,p1):
    newUp = om.MVector(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]).normal()
    quaternion = om.MQuaternion()
    if newUp.y > 0:
        quaternion = om.MVector.kYaxisVector.rotateTo(newUp)
    elif newUp.y < 0:
        quaternion = om.MVector.kYnegAxisVector.rotateTo(newUp)

    return quaternion


if __name__ == '__main__':
    bones_names = cmds.ls(selection=True)

    for bone_name in bones_names:
        selection = om.MSelectionList()
        selection.add(bone_name)
        bone_obj = selection.getDependNode(0)
        try:
            p0,p1,quaternion = get_joint_transform(bone_obj)
        except:
            raise
        else:
            def saveFile():
                dirPath = cmds.workspace(q=True, rootDirectory=True)
                with open(dirPath+'data/'+bone_name+'.dat','w') as f:
                    f.write('joint position: {} {} {}\n'.format(p0[0],p0[1],p0[2]))
                    f.write('next joint position: {} {} {}\n'.format(p1[0],p1[1],p1[2]))
                    f.write('quaternion: {} {} {} {}\n'.format(quaternion[0],quaternion[1],quaternion[2],quaternion[3]))
    
            saveFile()



