import maya.cmds as cmds

dirPath=cmds.workspace(q=True, rootDirectory=True )
dirPath+='scripts/'
cmds.loadPlugin(dirPath+"cross_section_extract_cmd.py")
grp=cmds.ls('source_cross_section_group',transforms=True)
if not grp:
    grp=cmds.group( em=True, name='source_cross_section_group' )
    

ankle = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for i in ankle:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightAnkle",md=40,mw=True)
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftAnkle",md=40,mw=True)
    cmds.parent( tmp, grp )
    

foot = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
for i in ankle:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightFoot",md=40,mw=True)
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftFoot",md=40,mw=True)
    cmds.parent( tmp, grp )
    
    
thigh=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
for i in thigh:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightThigh",md=40,mw=True)
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftThigh",md=40,mw=True)
    cmds.parent( tmp, grp )


knee = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for i in knee:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightKnee",md=40,mw=True) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftKnee",md=40,mw=True) 
    cmds.parent( tmp, grp )  
      
      
leg=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
#data_dir = cmds.workspace(q=True, rootDirectory=True )+'data/'
#myfile = open(data_dir+"LeftLeg_16_cross_section_curves.mel", "w")
#myfile.write("// cross-section curves of a leg\n")
#myfile.close()
for i in leg:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightLeg",md=40,mw=True) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftLeg",md=40,mw=True) 
    cmds.parent( tmp, grp )
    
    #input_file = '{0}Source_LeftLeg_cross_section_u_at_{1}_percentage_worldspace.dat'.format(data_dir,int(i*100))
    
    #f=open(input_file,'r') 
    #str = f.readlines()
    #f.close()
    #myfile = open(data_dir+"LeftLeg_16_cross_section_curves.mel", "a")
    #myfile.write('curve -d 3')
    #for s in str:
        #s=s.replace('\n','')
        #s='\n-p '+s       
        #myfile.write(s)
    #myfile.write(";\n")
    #myfile.close()

    
neck=[1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
for i in neck:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Neck",md=80,mw=True)
    cmds.parent( tmp, grp )


shoulder=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
for i in neck:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Shoulder",md=120,mw=True)
    cmds.parent( tmp, grp )


chest=[0.86,0.66,0.44,0.22,0]
for i in chest:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Chest",md=120,mw=True)
    cmds.parent( tmp, grp )


belly = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.25, 0.18, 0.1, 0]
for i in belly:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Belly",md=120,mw=True)
    cmds.parent( tmp, grp )

file_paths = [
'Source_Chest_cross_section_u_at_22_percentage.dat',
'Source_Chest_cross_section_u_at_66_percentage.dat',
'Source_Chest_cross_section_u_at_86_percentage.dat'
]
delete_points_list = {'Source_Chest_cross_section_u_at_22_percentage.dat':[34,37,83,86],
                      'Source_Chest_cross_section_u_at_66_percentage.dat':[30,39,81,90],
                      'Source_Chest_cross_section_u_at_86_percentage.dat':[28,40,80,92],
                     }
                                    
dir_path = cmds.workspace(fn=True)+'/data/'
for file_path in file_paths:
    vertices = []
    numPt = 0
    with open(dir_path + file_path, 'r') as f:
        content = f.readlines()
        for line in content:
            p = line.split()
            numPt += 1
            vertices.append(om.MVector(float(p[0]), float(p[1]), float(p[2])))
    # To rebuild the curve:
    # delete the unnecessary points
    # approximate the new list of points with trigonometric functions
    del vertices[delete_points_list.get(file_path)[2]+1:delete_points_list.get(file_path)[3]]
    del vertices[delete_points_list.get(file_path)[0]+1:delete_points_list.get(file_path)[1]]

    with open(dir_path + 'removed_' + file_path, "w") as f:
        for v in vertices:
            f.write('{} {} {}\n'.format(v[0],v[1],v[2]))
        # rebuild the curve for every cross-section
        center = curve_fitting.getCenter_2D(vertices)
        angles = curve_fitting.calculateAngle_2D(vertices, center)
        d_bar = curve_fitting.get_d_bar_2D(vertices, center)
        J, coe = curve_fitting.findJ_2D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, start_index=0)

        # add points to replace the deleted points
        delta_angle = 2 * math.pi / numPt
        for index in range(delete_points_list.get(file_path)[0], delete_points_list.get(file_path)[1]):
            angles.insert(index, delta_angle + angles[index])
        for index in range(delete_points_list.get(file_path)[2], delete_points_list.get(file_path)[3]):
            angles.insert(index, delta_angle + angles[index])
        curve_fitting.form_vertices_of_fragment_3D()

belly = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.25, 0.18, 0.1, 0]
for i in belly:
    tmp=cmds.crossSectionExtract(mu=i, mmn="source_male_mesh", mbn="Source_Belly", md=120, mw=True)
    cmds.parent( tmp, grp )


hip=[0.75,0.5,0.25,0.0]
for i in hip:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_Hip",md=120,mw=True)
    cmds.parent( tmp, grp )

    
arm=[0.35,0.4,0.5,0.55,0.6,0.7,0.8,0.9,0.95,1.0]#arm=[0.32,0.55,0.91,1]
for i in arm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightArm",md=40,mw=True) 
    cmds.parent( tmp, grp )

    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftArm",md=40,mw=True) 
    cmds.parent( tmp, grp )

foreArm=[0.0,0.1,0.2,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
for i in foreArm:
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_RightForeArm",md=40,mw=True) 
    cmds.parent( tmp, grp )
    
    tmp=cmds.crossSectionExtract(mu=i,mmn="source_male_mesh",mbn="Source_LeftForeArm",md=40,mw=True) 
    cmds.parent( tmp, grp )


