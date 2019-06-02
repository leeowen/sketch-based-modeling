# cross_section_cmd.py
#
# DESCRIPTION:
#
#   Produces the MEL "crossSectionExtract" command.
#
#   To use the command, select the joint and mesh that you want data for and
#   then type the command, being sure to use the -mu/-myUParameter flag to specify
#   the location on the bone where a ray will cast out and later intersect with the mesh 
#   to form the cross-section curve.
# 
#   For example:
#    
#     crossSectionExtract -mu 0.1 -mmn "source_male_mesh" -mbn "Source_LeftUpLeg"
#
#   The output is the cross section curve extracted from a mesh named "Source_male_meshShape", 
#   paralled with a bone named "Source_LeftUpLeg" at the location u=0.1.
#------------------------------------------------- 
import maya.api.OpenMaya as om
import maya.cmds as cmds
import numpy, math
import collections

def maya_useNewAPI():
    pass


def linear_interpolate_3D(p1,p2,t):
    p1=om.MVector(p1[0],p1[1],p1[2])
    p2=om.MVector(p2[0],p2[1],p2[2])
    p=t*p2+(1.0-t)*p1
    return p


def line_normal(p0,p1,p2):
    x0=p0[0]
    y0=p0[1]
    x1=p1[0]
    x2=p2[0]
    y1=p1[1]
    y2=p2[1]
    z=p0[2]
    slop=-(x1-x2)/(y1-y2)
    xp=x0-1.0
    yp=xp*slop+y0-slop*x0
    n=om.MVector(xp,yp,0)
    n.normalize()
    return n
 
    
Local_Frame_Tuple=collections.namedtuple("LocalFrame",['xAxis','yAxis','zAxis'])


mesh_name=''
mesh_dagPath = om.MDagPath()
bone_name=''  
nextBone_name=''
bone_dagPath = om.MDagPath()  
nextBone_dagPath = om.MDagPath() 
u_parameter=0.5    #default value

sList = om.MGlobal.getActiveSelectionList()
iter=om.MItSelectionList (sList, om.MFn.kMesh)
i=0
            
while( not iter.isDone()): 
    # RETRIEVE THE MESH
    mesh_dagPath = iter.getDagPath()
    meshFn=om.MFnMesh( mesh_dagPath)
    name=meshFn.fullPathName().split('|')
    mesh_name=name[-1]
    i+=1
    iter.next()
          
if i==0: 
    raise ValueError("No mesh or mesh transform specified!")
elif i>1:
    raise ValueError("Multiple meshes or mesh transforms specified!")
    
 
iter.reset()
iter=om.MItSelectionList (sList, om.MFn.kJoint)
i=0
while not iter.isDone(): 
    #RETRIEVE THE JOINT
    bone_dagPath = iter.getDagPath()
    boneFn=om.MFnDagNode( bone_dagPath)
    name=boneFn.fullPathName().split('|')
    bone_name=name[-1]
    i+=1
    iter.next()
               
if i==0: 
    raise ValueError("No bone or bone transform specified!")
elif i>1:
    raise ValueError("Multiple bones or bone transforms specified!")

               
# GET BONE DATA
if boneFn.childCount()>0:
    next_bone_obj=boneFn.child(0)
elif boneFn.parentCount():
    next_bone_obj=boneFn.parent(0)
            
nextBoneFn=om.MFnDagNode(next_bone_obj)
nextBone_name=nextBoneFn.fullPathName().split('|')
nextBone_name=nextBone_name[-1]
nextBone_dagPath = nextBoneFn.getPath()

bone_position=cmds.xform( bone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
bone_position=om.MVector(bone_position[0],bone_position[1],bone_position[2])
nextBone_position=cmds.xform( nextBone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
nextBone_position=om.MVector(nextBone_position[0],nextBone_position[1],nextBone_position[2])
ray_center=linear_interpolate_3D(bone_position,nextBone_position,u_parameter)        
    
# GET THE LOCAL FRAME ON THE CHOOSEN BONE
xAxis=(bone_position-nextBone_position).normalize()
yAxis=line_normal(ray_center,bone_position,nextBone_position)# this is not a real normal! just an intermiediate vector                                                                            
zAxis=xAxis^yAxis
zAxis.normalize()
yAxis=xAxis^zAxis
yAxis.normalize()
local_frame=Local_Frame_Tuple(xAxis,yAxis,zAxis)
       
#EXTRACT CROSS SECTION CURVES
vdiv=50 #the number of v division
meshFn=om.MFnMesh(mesh_dagPath)
raySource=om.MFloatPoint(ray_center)
print raySource
output=""
for i in range(0,vdiv):
    angle=2*math.pi*float(i)/float(vdiv)
    ray=yAxis.rotateBy(om.MQuaternion(angle,xAxis))
    rayDirection=om.MFloatVector(ray)

    try:
        hitPoint, hitRayParam, hitFace, hitTriangle, hitBary1, hitBary2 = meshFn.closestIntersection(raySource,rayDirection,om.MSpace.kWorld,9999,False,idsSorted=False,tolerance=0.001)    
    except:
        raise
    else:
        output+=str(hitPoint.x)+","+str(hitPoint.y)+","+str(hitPoint.z)+"\n"
    print output
               

     
        

        
