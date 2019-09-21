import maya.cmds as cmds

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")
grp=cmds.ls('source_cross_section_group',transforms=True)
if not grp:
    grp=cmds.group( em=True, name='source_cross_section_group' )

head=[0.9,0.62,0.22]
for i in head:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Head",md=80)
    cmds.parent( tmp, grp )

    
neck=[0.8,0.3]
for i in neck:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Neck",md=80)
    cmds.parent( tmp, grp )


chest=[0.96,0.86,0.66,0.22]
for i in chest:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Chest",md=120)
    cmds.parent( tmp, grp )

#cmds.crossSectionExtract(mu=0.42,mmn="source_male_mesh",mbn="Source_Chest")
#cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_Chest")
belly=[0.7,0.3,0.18]
for i in belly:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Belly",md=120)
    cmds.parent( tmp, grp )


hip=[0.75,0.5,0.25,0.0]
for i in hip:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Hip",md=120)
    cmds.parent( tmp, grp )


thigh=[0.2,0.3,0.6,0.83,0.91,1]
for i in thigh:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightThigh",md=40)
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftThigh",md=40)
    cmds.parent( tmp, grp )


leg=[0.05,0.1,0.45,0.52,0.94,1]
for i in leg:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightLeg",md=40) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftLeg",md=40) 
    cmds.parent( tmp, grp )

    
arm=[0.3,0.32,0.4,0.5,0.55,0.6,0.7,0.8,0.9,0.95,1.0]#arm=[0.32,0.55,0.91,1]
for i in arm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightArm",md=40) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftArm",md=40) 
    cmds.parent( tmp, grp )
    

foreArm=[0.1,0.2,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0]#foreArm=[0.1,0.35,0.5,0.6,0.9,1]
for i in foreArm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightForeArm",md=40) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftForeArm",md=40) 
    cmds.parent( tmp, grp )


