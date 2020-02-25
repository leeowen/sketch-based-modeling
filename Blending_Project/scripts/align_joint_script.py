# Before run this script, please make sure the joints are unparented already!
import maya.cmds as cmds


# To avoid issues like freezeTransform, recommend rotate pivot to attain the position
head_joint_position=cmds.xform("Source_Head",absolute=True,query=True,worldSpace=True,rotatePivot=True)
head_joint_position=om.MVector(head_joint_position[0],head_joint_position[1],head_joint_position[2])
hip_joint_position=cmds.xform("Source_Hip",absolute=True,query=True,worldSpace=True,rotatePivot=True)
hip_joint_position=om.MVector(hip_joint_position[0],hip_joint_position[1],hip_joint_position[2])
neck_joint_position=cmds.xform("Source_Neck",absolute=True,query=True,worldSpace=True,rotatePivot=True)
neck_joint_position=om.MVector(neck_joint_position[0],neck_joint_position[1],neck_joint_position[2])
chest_joint_position=cmds.xform("Source_Chest",absolute=True,query=True,worldSpace=True,rotatePivot=True)
chest_joint_position=om.MVector(chest_joint_position[0],chest_joint_position[1],chest_joint_position[2])
top_of_head_joint_position=cmds.xform("Source_Top_Of_Head",absolute=True,query=True,worldSpace=True,rotatePivot=True)
top_of_head_joint_position=om.MVector(top_of_head_joint_position[0],top_of_head_joint_position[1],top_of_head_joint_position[2])
belly_joint_position=cmds.xform("Source_Belly",absolute=True,query=True,worldSpace=True,rotatePivot=True)
belly_joint_position=om.MVector(belly_joint_position[0],belly_joint_position[1],belly_joint_position[2])

line = top_of_head_joint_position - hip_joint_position

def new_joint_position(joint_position):
    t = (joint_position[1] - hip_joint_position[1])/line[1]
    return t
    
def linear_interpolate_3D(t):
    px=(1.0-t)*hip_joint_position[0]+t*top_of_head_joint_position[0]
    py=(1.0-t)*hip_joint_position[1]+t*top_of_head_joint_position[1]
    pz=(1.0-t)*hip_joint_position[2]+t*top_of_head_joint_position[2]
    p=om.MVector(px,py,pz)
    return p

neck_joint_position = linear_interpolate_3D(new_joint_position(neck_joint_position))
cmds.setAttr("Source_Neck"+'.translate',neck_joint_position[0],neck_joint_position[1],neck_joint_position[2])
head_joint_position = linear_interpolate_3D(new_joint_position(head_joint_position))
cmds.setAttr("Source_Head"+'.translate',head_joint_position[0],head_joint_position[1],head_joint_position[2])
chest_joint_position = linear_interpolate_3D(new_joint_position(chest_joint_position))
cmds.setAttr("Source_Chest"+'.translate',chest_joint_position[0],chest_joint_position[1],chest_joint_position[2])
belly_joint_position = linear_interpolate_3D(new_joint_position(belly_joint_position))
cmds.setAttr("Source_Belly"+'.translate',belly_joint_position[0],belly_joint_position[1],belly_joint_position[2])
