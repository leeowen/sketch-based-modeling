import maya.api.OpenMaya as om
import os
sys.path.append(cmds.workspace(fn=True)+'/scripts/')
import curve_fitting
import math


file_paths = [
'modified_Source_Chest_cross_section_u_at_0_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_22_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_66_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_86_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_100_percentage_worldspace.dat'
]

Ea_criteria = 0.01
Em_criteria = 0.025
dirPath = cmds.workspace(fn=True)+'/data/'
for file_path in file_paths:
    vertices = []
    numPt = 0
    Ea = 0
    Em = 0
    with open(dirPath + file_path, 'r') as f:
        content = f.readlines()
        for line in content:
            p = line.split()
            numPt += 1
            vertices.append(om.MVector(float(p[0]), float(p[1]), float(p[2])))
    center = curve_fitting.getCenter_3D(vertices)
    file_name = file_path.replace('Source_', '')
    directory = dirPath + 'one_piece_' + file_name.split('_cross_section_')[0]

    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    save_file_path = directory + '/' + file_name

    isClosed = True

    #curve_fitting.maya_polygon_plane(vertices)

    center = curve_fitting.getCenter_3D(vertices)
    angles = [i * 2 * math.pi / numPt for i in range(numPt)]
    #angles = curve_fitting.calculateAngle_3D(vertices, center)
    d_bar = curve_fitting.get_d_bar_3D(vertices, center)

    J, coe = curve_fitting.findJ_3D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, 0)
    new_vertices = curve_fitting.form_vertices_of_fragment_3D(coe, center, angles)
    Ea, Em = curve_fitting.calculate_Ea_Em_3D(vertices, new_vertices, d_bar, 0)
    print file_name, J, Ea, Em
    with open(save_file_path, "w+") as f:
        f.write('range:0-360\n')
        f.write('center: {} {} {}\n'.format(center[0], center[1], center[2]))
        f.write('a: ')
        for a in coe[0]:
            f.write(str(a) + ' ')
        f.write('\n')
        f.write('b: ')
        for b in coe[1]:
            f.write(str(b) + ' ')
        f.write('\n')
        f.write('c: ')
        for c in coe[2]:
            f.write(str(c) + ' ')
        f.write('\n')
        f.write('angles: ')
        for angle in angles:
            f.write(str(angle) + ' ')
        f.write('\n')
