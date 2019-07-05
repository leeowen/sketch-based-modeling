import maya.cmds as cmds

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")

grp=cmds.group( em=True, name='source_cross_section_group' )
grp_meta=cmds.group( em=True, name='source_meta_cross_section_group' )

head=[0.9,0.62,0.22]
for i in head:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Head",md=40)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
neck=[0.8,0.3]
for i in neck:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Neck",md=40)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

chest=[0.96,0.86,0.66,0.22]
for i in chest:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Chest",md=80)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
#cmds.crossSectionExtract(mu=0.42,mmn="source_male_mesh",mbn="Source_Chest")
#cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_Chest")
belly=[0.7,0.3,0.18]
for i in belly:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Belly",md=80)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

hip=[0.75,0.5,0.25,0.0]
for i in hip:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Hip",md=80)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

thigh=[0.2,0.3,0.6,0.83,0.91,1]
for i in thigh:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightThigh",md=40)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftThigh",md=40)
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

leg=[0.05,0.1,0.45,0.52,0.94,1]
for i in leg:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightLeg",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftLeg",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
arm=[0.32,0.55,0.91,1]
for i in arm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightArm",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftArm",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

foreArm=[0.1,0.35,0.5,0.6,0.9,1]
for i in foreArm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightForeArm",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftForeArm",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)

grp=cmds.group( em=True, name='target_cross_section_group' )
arm=[0.32,0.55,0.91,1]
for i in arm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="target_male_mesh",mbn="Target_RightArm",md=40) 
    cmds.parent( tmp[0], grp )
    cmds.parent(tmp[1],grp_meta)