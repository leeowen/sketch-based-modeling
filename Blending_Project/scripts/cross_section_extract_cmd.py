# cross_section_extract_cmd.py
#
# DESCRIPTION:
#
#   Produces the MEL "crossSectionExtract" command.
#
#   To use the command, select the joint and mesh that you want data for and
#   then type the command, being sure to use the -mu/-myUParameter flag to specify
#   the location on the Joint where a ray will cast out and later intersect with the mesh 
#   to form the cross-section curve.
# 
#   For example:
#    
#     crossSectionExtract -mu 0.1 -md 40 -mmn "source_male_mesh" -mbn "Source_LeftUpLeg"
#
#   The output is made up by two parts. One is a cross section curve extracted from mesh "Source_male_meshShape", 
#   paralled with a Joint named "Source_LeftUpLeg" at the location u=0.1, and it is in world space. 
#   The other is a file contains the extracted cross section curve's points' position informations, which are in object space.
#------------------------------------------------- 
import maya.api.OpenMaya as om
import maya.cmds as cmds
import math
import os.path

def maya_useNewAPI():
    pass

##############################################################################
##
## Command class implementation
##
##############################################################################
class CrossSectionExtractCmd(om.MPxCommand):
    kPluginName='crossSectionExtract'
    def __init__(self):
        om.MPxCommand.__init__(self)
        self.mesh_name=''
        self.mesh_dagPath = om.MDagPath()
        self.joint_name=''  
        self.nextJoint_name=''
        self.Joint_dagPath = om.MDagPath()  
        self.nextJoint_dagPath = om.MDagPath() 
        self.u_parameter=0.0    #default value
        self.division=20

        
    @staticmethod
    def creator():
        return CrossSectionExtractCmd()
        
        
    @staticmethod
    def newSyntax():
        syntax=om.MSyntax()
        syntax.addFlag('-mu','-myUparameter',om.MSyntax.kDouble)
        syntax.addFlag('-mmn','-myMeshName',om.MSyntax.kString)
        syntax.addFlag('-mbn','-myJointName',om.MSyntax.kString)
        syntax.addFlag('-md','-myDivision',om.MSyntax.kLong)
        return syntax
        
    
    def isUndoable(self):
        ''' This function indicates whether or not the command is undoable. If the
        command is undoable, each executed instance of that command is added to the
        undo queue. '''
        
        # We must return True to specify that this command is undoable.
        return False

                             
    def parseArguments(self, args):
        argData=om.MArgDatabase(self.syntax(),args)
        if argData.isFlagSet('-mu'):
            self.u_parameter = argData.flagArgumentDouble( '-mu', 0 )
        if argData.isFlagSet('-mmn'):
            self.mesh_name = argData.flagArgumentString( '-mmn', 0 )
        if argData.isFlagSet('-mbn'):
            self.joint_name = argData.flagArgumentString( '-mbn', 0 )  
        if argData.isFlagSet('-md'):
            self.division = argData.flagArgumentInt( '-md', 0 )  
                
        # WHEN NO MESH IS SPECIFIED IN THE COMMAND, GET THE FIRST SELECTED MESH FROM THE SELECTION LIST:
        if (self.mesh_name == "") or (self.joint_name == ""):
            sList = om.MGlobal.getActiveSelectionList()
                      
            iter=om.MItSelectionList (sList, om.MFn.kMesh)
            i=0
            
            while( not iter.isDone()): 
                # RETRIEVE THE MESH
                self.mesh_dagPath = iter.getDagPath()
                meshFn=om.MFnMesh(self.mesh_dagPath)
                name=meshFn.fullPathName().split('|')
                self.mesh_name=name[-1]
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
                self.Joint_dagPath = iter.getDagPath()
                jointFn=om.MFnDagNode(self.Joint_dagPath)
                name=jointFn.name()
                i+=1
                iter.next()
                
            if i==0: 
                raise ValueError("No Joint or Joint transform specified!")
            elif i>1:
                raise ValueError("Multiple Joints or Joint transforms specified!")

         
        # get dagpath from mesh's and joint's names 
        else:   
            try:
                selectionList = om.MGlobal.getSelectionListByName(self.mesh_name)
            except:
                raise
            else:
                self.mesh_dagPath = selectionList.getDagPath( 0 )
            try:
                selectionList = om.MGlobal.getSelectionListByName(self.joint_name)
            except:
                raise
            else:
                self.Joint_dagPath = selectionList.getDagPath( 0 )


    def doIt(self,args):
        """
        The doIt method should collect whatever information is required to do the task, and store it in local class data. 
        It should finally call redoIt to make the command happen. 
        """
        self.parseArguments(args)
        
        # GET Joint DATA
        jointFn=om.MFnDagNode(self.Joint_dagPath)
        if jointFn.childCount()==1:
            next_joint_obj=jointFn.child(0)
        elif jointFn.childCount()>1:# Belly, Neck
            for i in range(jointFn.childCount()):
                tmpObj=jointFn.child(i)
                tmpFn=om.MFnDagNode(tmpObj)
                tmpName=tmpFn.name()
                if ("Chest" in tmpName) or ("Head" in tmpName):
                    next_joint_obj=jointFn.child(i)
        elif jointFn.childCount()==0 and jointFn.parentCount()==1:# Hip Hand Toe Top
            next_joint_obj=jointFn.parent(0)
            
        nextJointFn=om.MFnDagNode(next_joint_obj)
        self.nextJoint_dagPath = nextJointFn.getPath()
        self.nextJoint_name=nextJointFn.name()
                       
        # To avoid issues like freezeTransform, recommend rotate pivot to attain the position
        joint_position=cmds.xform(self.joint_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        joint_position=om.MVector(joint_position[0],joint_position[1],joint_position[2])
        nextjoint_position=cmds.xform(self.nextJoint_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        nextjoint_position=om.MVector(nextjoint_position[0],nextjoint_position[1],nextjoint_position[2])
        ray_center=linear_interpolate_3D(joint_position,nextjoint_position,self.u_parameter)                
       
       # GET THE LOCAL FRAME ON THE CHOOSEN JOINT
        VAxis=(joint_position-nextjoint_position).normalize()
        if VAxis.y<0:
            VAxis=-VAxis                                                                          
        UAxis=om.MVector.kXaxisVector # this is not a real normal! just an intermiediate vector 
        if abs(VAxis*UAxis)>0.99:
            UAxix=om.MVector.kZaxisVector
        WAxis=UAxis^VAxis
        WAxis.normalize()
        UAxis=VAxis^WAxis
        UAxis.normalize()

        # Find the quaternion that form the coordinate transformation of 
        # joint's local framework: how the new coordinate get back into world coordinate
        quaternion=getQuaternion(UAxis,VAxis,WAxis)
        
        transform_name=self.joint_name+"_cross_section"+"_u_at_"+str(int(self.u_parameter*100))+"_percentage"
        dirPath=cmds.workspace(q=True, rootDirectory=True )
               
        #EXTRACT CROSS SECTION CURVES
        meshFn=om.MFnMesh(self.mesh_dagPath)
        raySource=om.MFloatPoint(ray_center)
        eps=om.MPointArray()
        eps.clear()
        epsLocal=om.MPointArray()
        epsLocal.clear()
            
        for i in range(0,self.division):
            angle=2*math.pi*float(i)/float(self.division)
            ray=UAxis.rotateBy(om.MQuaternion(angle,VAxis))
            rayDirection=om.MFloatVector(ray)
            
            try:
                hitPoint, hitRayParam, hitFace, hitTriangle, hitBary1, hitBary2 = meshFn.closestIntersection(raySource,rayDirection,om.MSpace.kWorld,9999,False,idsSorted=False,tolerance=0.001)    
            except:
                raise
            else:
                eps.append(hitPoint)
                hitPointLocal=om.MVector(hitPoint[0],hitPoint[1],hitPoint[2])
                hitPointLocal-=ray_center
                # This is the quaternion that transform a point from world coordinate into local coordinate 
                # Hence, we need to inverse the quaternion we get above
                hitPointLocal=hitPointLocal.rotateBy(quaternion.inverse())
                epsLocal.append(hitPointLocal)
                
        curveFn=om.MFnNurbsCurve()
        self.cross_section_obj=curveFn.createWithEditPoints(eps,2,om.MFnNurbsCurve.kClosed, False, False, True)
       
        dgFn=om.MFnDependencyNode(self.cross_section_obj)
        dgFn.setName(transform_name)
                
        # add custom attribute to node
        attrFn=om.MFnNumericAttribute()
        uAttr=attrFn.create("uParameter","u",om.MFnNumericData.kFloat,self.u_parameter)
        attrFn.readable=True 
        attrFn.storable=True # fairly consistent, won't change in compute() or get updated by upstream node, etc
        attrFn.writable=True 
        attrFn.keyable=True
        attrFn.hidden=False
        dgFn.addAttribute( uAttr)
                       
        om.MPxCommand.setResult(om.MFnDependencyNode(self.cross_section_obj).name())
                                #,om.MFnDependencyNode(self.cross_section_obj_meta).name()])
        save_path=dirPath+'data/'
        complete_name=os.path.join(save_path,transform_name+'.dat')
        outFile=open(complete_name,'w+')# Here we used "w" letter in our argument, which indicates write and the plus sign that means it will create a file if it does not exist in library
        for i in range(0,self.division):
            x=epsLocal[i].x
            y=epsLocal[i].y
            z=epsLocal[i].z
            outFile.write('{} {} {} \n'.format(x,y,z))
            
        outFile.close()
        
def initializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.registerCommand(CrossSectionExtractCmd.kPluginName, CrossSectionExtractCmd.creator, CrossSectionExtractCmd.newSyntax)
    except:
        raise
        
        
def uninitializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.deregisterCommand(CrossSectionExtractCmd.kPluginName)
    except:
        raise
        
        
def linear_interpolate_3D(p1,p2,t):
    p1=om.MVector(p1[0],p1[1],p1[2])
    p2=om.MVector(p2[0],p2[1],p2[2])
    p=t*p2+(1.0-t)*p1
    return p
    

def getQuaternion(UAxis,VAxis,WAxis):
    q=om.MQuaternion()
    qy=om.MVector.kYaxisVector.rotateTo(VAxis)
    q=qy
    xRotated=om.MVector.kXaxisVector.rotateBy(q)
    angle=math.acos(xRotated*UAxis)
    qx=om.MQuaternion(angle,VAxis)
    if not UAxis.isEquivalent(xRotated.rotateBy(qx),1.0e-5):
        angle=2*math.pi-angle
        qx=om.MQuaternion(angle,VAxis)
        
    q=q*qx
    return q


