# get all resources from Maya

import maya.cmds as cmds

for item in cmds.resourceManager(nameFilter='*png'):
    print(item)