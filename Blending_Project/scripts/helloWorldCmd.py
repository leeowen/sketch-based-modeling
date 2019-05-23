import maya.api.OpenMaya as om

def maya_useNewAPI():
    pass


class helloWorldCmd(om.MPxCommand):
    kPluginName='helloWorldCmd'
    def __init__(self):
        om.MPxCommand.__init__(self)
        
    @staticmethod
    def creator():
        return helloWorldCmd()
        
    @staticmethod
    def syntaxCreator():
        syntax=om.MSyntax()
        syntax.addFlag('-mu','-myUparameter',om.MSyntax.kDouble)
        return syntax
        
    def doIt(self,args):
        self.parseArguments(args)
        print("hello world with syntax -u:",self.u)


    def parseArguments(self, args):
        argData=om.MArgDatabase(self.syntax(),args)
        if argData.isFlagSet('-mu'):
            self.u = argData.flagArgumentDouble( '-mu', 0 )

def initializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.registerCommand(helloWorldCmd.kPluginName, helloWorldCmd.creator, helloWorldCmd.syntaxCreator)
    except:
        raise
        
        
def uninitializePlugin(plugin):
    pluginFn = om.MFnPlugin(plugin)
    try:
        pluginFn.deregisterCommand(helloWorldCmd.kPluginName)
    except:
        raise