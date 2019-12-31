import maya.cmds as cmds

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")
grp=cmds.ls('source_cross_section_group',transforms=True)
if not grp:
    grp=cmds.group( em=True, name='source_cross_section_group' )
    

#thigh=[0.2,0.3,0.45,0.6,0.7,0.83,0.91,1]
#for i in thigh:
    #tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightThigh",md=40,mw=True)
    #cmds.parent( tmp, grp )
    
    #tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftThigh",md=40,mw=True)
    #cmds.parent( tmp, grp )


leg=[0.1,0.15,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,0.95,1.0]
data_dir = cmds.workspace(q=True, rootDirectory=True )+'data/'
myfile = open(data_dir+"LeftLeg_16_cross_section_curves.mel", "w")
myfile.write("// cross-section curves of a leg\n")
myfile.close()
for i in leg:
    #tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightLeg",md=40,mw=True) 
    #cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftLeg",md=40,mw=True) 
    cmds.parent( tmp, grp )
    
    input_file = '{0}Source_LeftLeg_cross_section_u_at_{1}_percentage_worldspace.dat'.format(data_dir,int(i*100))
    
    f=open(input_file,'r') 
    str = f.readlines()
    f.close()
    myfile = open(data_dir+"LeftLeg_16_cross_section_curves.mel", "a")
    myfile.write('curve -d 3')
    for s in str:
        s=s.replace('\n','')
        s='\n-p '+s       
        myfile.write(s)
    myfile.write(";\n")
    myfile.close()
        
        
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

    
arm=[0.4,0.5,0.55,0.6,0.7,0.8,0.9,0.95,1.0]#arm=[0.32,0.55,0.91,1]
for i in arm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightArm",md=40) 
    cmds.parent( tmp, grp )

leftArm=[0.16,0.2,0.3,0.4,0.50,0.6,0.7,0.8,0.9,0.95,1.0]
for i in leftArm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftArm",md=40) 
    cmds.parent( tmp, grp )

foreArm=[0.0,0.1,0.2,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0]#foreArm=[0.1,0.35,0.5,0.6,0.9,1]
for i in foreArm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightForeArm",md=40) 
    cmds.parent( tmp, grp )
    
    #tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftForeArm",md=40) 
    #cmds.parent( tmp, grp )


