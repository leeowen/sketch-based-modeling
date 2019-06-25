import maya.cmds as cmds

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")

cmds.crossSectionExtract(mu=0.9,mmn="source_male_mesh",mbn="Source_Head")
cmds.crossSectionExtract(mu=0.62,mmn="source_male_mesh",mbn="Source_Head")
cmds.crossSectionExtract(mu=0.22,mmn="source_male_mesh",mbn="Source_Head")

cmds.crossSectionExtract(mu=0.8,mmn="source_male_mesh",mbn="Source_Neck")
cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_Neck")

cmds.crossSectionExtract(mu=0.96,mmn="source_male_mesh",mbn="Source_Chest")
cmds.crossSectionExtract(mu=0.86,mmn="source_male_mesh",mbn="Source_Chest")
cmds.crossSectionExtract(mu=0.66,mmn="source_male_mesh",mbn="Source_Chest")
cmds.crossSectionExtract(mu=0.42,mmn="source_male_mesh",mbn="Source_Chest")
cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_Chest")
cmds.crossSectionExtract(mu=0.22,mmn="source_male_mesh",mbn="Source_Chest")

cmds.crossSectionExtract(mu=0.7,mmn="source_male_mesh",mbn="Source_Belly")
cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_Belly")
cmds.crossSectionExtract(mu=0.18,mmn="source_male_mesh",mbn="Source_Belly")

cmds.crossSectionExtract(mu=0.75,mmn="source_male_mesh",mbn="Source_Hip")
cmds.crossSectionExtract(mu=0.5,mmn="source_male_mesh",mbn="Source_Hip")
cmds.crossSectionExtract(mu=0.25,mmn="source_male_mesh",mbn="Source_Hip")
cmds.crossSectionExtract(mu=0.0,mmn="source_male_mesh",mbn="Source_Hip")

cmds.crossSectionExtract(mu=0.2,mmn="source_male_mesh",mbn="Source_RightUpLeg")
cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_RightUpLeg")
cmds.crossSectionExtract(mu=0.6,mmn="source_male_mesh",mbn="Source_RightUpLeg")
cmds.crossSectionExtract(mu=0.83,mmn="source_male_mesh",mbn="Source_RightUpLeg") 
cmds.crossSectionExtract(mu=0.91,mmn="source_male_mesh",mbn="Source_RightUpLeg") 
cmds.crossSectionExtract(mu=1,mmn="source_male_mesh",mbn="Source_RightUpLeg") 

cmds.crossSectionExtract(mu=0.05,mmn="source_male_mesh",mbn="Source_RightLeg") 
cmds.crossSectionExtract(mu=0.1,mmn="source_male_mesh",mbn="Source_RightLeg") 
cmds.crossSectionExtract(mu=0.45,mmn="source_male_mesh",mbn="Source_RightLeg")
cmds.crossSectionExtract(mu=0.52,mmn="source_male_mesh",mbn="Source_RightLeg")
cmds.crossSectionExtract(mu=0.94,mmn="source_male_mesh",mbn="Source_RightLeg")
cmds.crossSectionExtract(mu=1,mmn="source_male_mesh",mbn="Source_RightLeg")

cmds.crossSectionExtract(mu=0.2,mmn="source_male_mesh",mbn="Source_LeftUpLeg")
cmds.crossSectionExtract(mu=0.3,mmn="source_male_mesh",mbn="Source_LeftUpLeg")
cmds.crossSectionExtract(mu=0.6,mmn="source_male_mesh",mbn="Source_LeftUpLeg")
cmds.crossSectionExtract(mu=0.83,mmn="source_male_mesh",mbn="Source_LeftUpLeg") 
cmds.crossSectionExtract(mu=0.91,mmn="source_male_mesh",mbn="Source_LeftUpLeg") 
cmds.crossSectionExtract(mu=1,mmn="source_male_mesh",mbn="Source_LeftUpLeg") 

cmds.crossSectionExtract(mu=0.05,mmn="source_male_mesh",mbn="Source_LeftLeg") 
cmds.crossSectionExtract(mu=0.1,mmn="source_male_mesh",mbn="Source_LeftLeg") 
cmds.crossSectionExtract(mu=0.45,mmn="source_male_mesh",mbn="Source_LeftLeg")
cmds.crossSectionExtract(mu=0.52,mmn="source_male_mesh",mbn="Source_LeftLeg")
cmds.crossSectionExtract(mu=0.94,mmn="source_male_mesh",mbn="Source_LeftLeg")
cmds.crossSectionExtract(mu=1,mmn="source_male_mesh",mbn="Source_LeftLeg")