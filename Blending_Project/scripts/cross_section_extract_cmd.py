# cross_section_cmd.py
# Author: Xiaosong Yang, Ouwen Li
# Email: leeowen988@gmail.com
# Version: 2.0
# Autodesk Maya Command template API 2.0
# Copyright (C) 2019 Xiaosong Yang
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
        self.bone_name=''    
        self.u_parameter=0.0    #default value
        
    @staticmethod
    def creator():
        return CrossSectionExtractCmd()
        
    @staticmethod
    def syntaxCreator():
        syntax=om.MSyntax()
        syntax.addFlag('-mu','-myUparameter',om.MSyntax.kDouble)
        syntax.addFlag('-mmn','-myMeshName',om.MSyntax.kDouble)
        syntax.addFlag('-mbn','-myBoneName',om.MSyntax.kDouble)
        return syntax
        
    def doIt(self,args):
        self.parseArguments(args)
        print("cross section extract at u=",self.u)


    def parseArguments(self, args):
        argData=om.MArgDatabase(self.syntax(),args)
        if argData.isFlagSet('-mu'):
            self.u = argData.flagArgumentDouble( '-mu', 0 )
        elif argData.isFlagSet('-mmn'):
            self.mesh_name = argData.flagArgumentDouble( '-mmn', 0 )
            

def initializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.registerCommand(CrossSectionExtractCmd.kPluginName, CrossSectionExtractCmd.creator, CrossSectionExtractCmd.syntaxCreator)
    except:
        raise
        
        
def uninitializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.deregisterCommand(CrossSectionExtractCmd.kPluginName)
    except:
        raise



        
