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

##############################################################################
##
## Command class implementation
##
##############################################################################

Local_Frame_Tuple=collections.namedtuple("LocalFrame",['xAxis','yAxis','zAxis'])

class CrossSectionExtractCmd(om.MPxCommand):
    kPluginName='crossSectionExtract'
    def __init__(self):
        om.MPxCommand.__init__(self)
        self.mesh_name=''
        self.mesh_dagPath = om.MDagPath()
        self.bone_name=''  
        self.nextBone_name=''
        self.bone_dagPath = om.MDagPath()  
        self.nextBone_dagPath = om.MDagPath() 
        self.u_parameter=0.0    #default value

        
    @staticmethod
    def creator():
        return CrossSectionExtractCmd()
        
        
    @staticmethod
    def newSyntax():
        syntax=om.MSyntax()
        syntax.addFlag('-mu','-myUparameter',om.MSyntax.kDouble)
        syntax.addFlag('-mmn','-myMeshName',om.MSyntax.kString)
        syntax.addFlag('-mbn','-myBoneName',om.MSyntax.kString)
        return syntax
        
    
    def isUndoable(self):
        ''' This function indicates whether or not the command is undoable. If the
        command is undoable, each executed instance of that command is added to the
        undo queue. '''
        
        # We must return True to specify that this command is undoable.
        return True
        
            
    def redoIt(self):   
        bone_position=cmds.xform(self.bone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        bone_position=om.MVector(bone_position[0],bone_position[1],bone_position[2])
        nextBone_position=cmds.xform(self.nextBone_name,absolute=True,query=True,worldSpace=True,rotatePivot=True)
        nextBone_position=om.MVector(nextBone_position[0],nextBone_position[1],nextBone_position[2])
        ray_center=linear_interpolate_3D(bone_position,nextBone_position,self.u_parameter)        
        

        # GET THE LOCAL FRAME ON THE CHOOSEN BONE
        xAxis=(bone_position-nextBone_position).normalize()
        yAxis=line_normal(ray_center,bone_position,nextBone_position)# this is not a real normal! just an intermiediate vector                                                                            
        zAxis=xAxis^yAxis
        zAxis.normalize()
        yAxis=xAxis^zAxis
        yAxis.normalize()
        local_frame=Local_Frame_Tuple(xAxis,yAxis,zAxis)
        
        #EXTRACT CROSS SECTION CURVES
        vdiv=20 #the number of v division
        meshFn=om.MFnMesh(self.mesh_dagPath)
        raySource=om.MFloatPoint(ray_center)
        #print raySource
        eps=om.MPointArray()
        eps.clear()
        for i in range(0,vdiv):
            angle=2*math.pi*float(i)/float(vdiv)
            ray=yAxis.rotateBy(om.MQuaternion(angle,xAxis))
            rayDirection=om.MFloatVector(ray)
        
            try:
                hitPoint, hitRayParam, hitFace, hitTriangle, hitBary1, hitBary2 = meshFn.closestIntersection(raySource,rayDirection,om.MSpace.kWorld,9999,False,idsSorted=False,tolerance=0.001)    
            except:
                raise
            else:
                eps.append(hitPoint)
                
        curveFn=om.MFnNurbsCurve()
        self.cross_section_obj=curveFn.createWithEditPoints(eps,2,om.MFnNurbsCurve.kClosed, False, False, True)
                 
        transform_name=self.bone_name+"_cross_section"
        curve_name="u_at_"+str(int(self.u_parameter*100))+"_percentage"
        
        dgFn=om.MFnDependencyNode(self.cross_section_obj)
        dgFn.setName(transform_name+'_'+curve_name)
              

    def undoIt(self):
        dagFn=om.MFnDagNode(self.cross_section_obj)
        child_curve=dagFn.child(0)

        dgModifier=om.MDGModifier()
        
        # always delete the child before deleting parent node
        dgModifier.deleteNode(child_curve)
        
        dgModifier.deleteNode(self.cross_section_obj)

        dgModifier.doIt()
        
           
    def parseArguments(self, args):
        argData=om.MArgDatabase(self.syntax(),args)
        if argData.isFlagSet('-mu'):
            self.u_parameter = argData.flagArgumentDouble( '-mu', 0 )
        if argData.isFlagSet('-mmn'):
            self.mesh_name = argData.flagArgumentString( '-mmn', 0 )
        if argData.isFlagSet('-mbn'):
            self.bone_name = argData.flagArgumentString( '-mbn', 0 )  
            
        # WHEN NO MESH IS SPECIFIED IN THE COMMAND, GET THE FIRST SELECTED MESH FROM THE SELECTION LIST:
        if (self.mesh_name == "") or (self.bone_name == ""):
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
                self.bone_dagPath = iter.getDagPath()
                boneFn=om.MFnDagNode(self.bone_dagPath)
                name=boneFn.fullPathName().split('|')
                self.bone_name=name[-1]
                i+=1
                iter.next()
                
            if i==0: 
                raise ValueError("No bone or bone transform specified!")
            elif i>1:
                raise ValueError("Multiple bones or bone transforms specified!")

         
        # get dagpath from mesh's and joint's names 
        else:   
            try:
                selectionList = om.MGlobal.getSelectionListByName(self.mesh_name)
            except:
                raise
            else:
                self.mesh_dagPath = selectionList.getDagPath( 0 )
            try:
                selectionList = om.MGlobal.getSelectionListByName(self.bone_name)
            except:
                raise
            else:
                self.bone_dagPath = selectionList.getDagPath( 0 )


    def doIt(self,args):
        """
        The doIt method should collect whatever information is required to do the task, and store it in local class data. 
        It should finally call redoIt to make the command happen. 
        """
        self.parseArguments(args)
        
        # GET BONE DATA
        boneFn=om.MFnDagNode(self.bone_dagPath)
        if boneFn.childCount()>0:
            next_bone_obj=boneFn.child(0)
        elif boneFn.parentCount():
            next_bone_obj=boneFn.parent(0)
            
        nextBoneFn=om.MFnDagNode(next_bone_obj)
        self.nextBone_name=nextBoneFn.fullPathName().split('|')
        self.nextBone_name=self.nextBone_name[-1]
        self.nextBone_dagPath = nextBoneFn.getPath()
        
        self.redoIt()                                            


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
        
