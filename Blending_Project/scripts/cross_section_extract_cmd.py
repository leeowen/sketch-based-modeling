# cross_section_cmd.py
# Author: Xiaosong Yang, Ouwen Li
# Email: leeowen988@gmail.com
# Version: 2.0
# Autodesk Maya Command template API 2.0
# Copyright (C) 2019 Xiaosong Yang
#-------------------------------------------------
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
#     crossSectionExtract -mu 0.1 -mmn "Source_male_meshShape" -mbn "Source_LeftUpLeg"
#
#   The output is the cross section curve extracted from a mesh named "Source_male_meshShape", 
#   paralled with a bone named "Source_LeftUpLeg" at the location u=0.1.
#------------------------------------------------- 
import maya.api.OpenMaya as om

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
        self.bone_name=''  
        self.bone_dagPath = om.MDagPath()  
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
        
        
    def doIt(self,args):
        """
        The doIt method should collect whatever information is required to do the task, and store it in local class data. 
        It should finally call redoIt to make the command happen. 
        """
        self.parseArguments(args)
        
        # PREP THE DATA
        boneFn=om.MFnDagNode(self.bone_dagPath)
        next_bone_obj=om.MObject()
        selList=om.MSelectionList()
        selList.add(self.bone_dagPath)
        selList.getDependNode(0,next_bone_obj)
        if boneFn.hasChild():
            next_bone_obj=boneFn.child(0)
        elif boneFn.hasParent():
            next_bone_obj=boneFn.parent(0)
        print next_bone_obj
                
        #DO THE WORK
        self.redoIt()


    def parseArguments(self, args):
        argData=om.MArgDatabase(self.syntax(),args)
        if argData.isFlagSet('-mu'):
            self.u = argData.flagArgumentDouble( '-mu', 0 )
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
                i+=1
                iter.next()
                
            if i==0: 
                raise ValueError("No bone or bone transform specified!")
            elif i>1:
                raise ValueError("Multiple bones or bone transforms specified!")

         
        # get dagpath from mesh's and joint's names 
        else:   
            print self.bone_dagPath
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

                                      
        
    def redoIt(self): 
        """
        The redoIt method should do the actual work, using only the local class data. 
        """
        pass    
        
        
    def undoIt(self):
        """
        The undoIt method should undo the actual work, again using only the local class data.
        """
        pass                 
                               

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



        
