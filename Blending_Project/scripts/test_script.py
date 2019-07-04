import maya.api.OpenMaya as om
import maya.cmds as cmds
import math,sys
sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy
def maya_useNewAPI():
    pass


quaternion=om.MQuaternion()
print quaternion

cy=cmds.cylinder()
cyName=cy[0]
sList=om.MSelectionList()
sList.add(cyName)
cyDag=sList.getDagPath(0)
cyObj=sList.getDependNode(0)

quaternion0=om.MQuaternion(1,1,1,math.cos(math.radians(60)))
transformFn=om.MFnTransform(cyObj)
transformFn.rotateBy(quaternion0,om.MSpace.kObject)
translate0=om.MVector(2,2,2)
print translate0
transformFn.translateBy(translate0,om.MSpace.kWorld)
quaternion1=om.MQuaternion(0,1,0,math.cos(math.radians(-60)))
transformFn.rotateBy(quaternion1,om.MSpace.kObject)

dgFn=om.MFnDagNode(cyObj)
# Add MPointArray attribute
tAttrFn=om.MFnTypedAttribute()
tAttr=tAttrFn.create('points','p',om.MFnData.kPointArray)
tAttrFn.array=True
tAttrFn.readable=True
tAttrFn.writable=True
tAttrFn.storable=True
dgFn.addAttribute(tAttr)

arrayPlug=dgFn.findPlug('points',False)

       
