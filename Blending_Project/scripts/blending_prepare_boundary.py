import maya.api.OpenMaya as om
import os
sys.path.append(cmds.workspace(fn=True)+'/scripts/')
import curve_fitting


dirPath = cmds.workspace(fn=True)+'/data/'
#source_file_path = 'Belly/Belly_cross_section_u_at_100_percentage_worldspace.dat'
#save_file_path = dirPath + 'Chest/Chest_cross_section_u_at_0_percentage_worldspace.dat'
#cut_points = [27, 42, 78, 93]
source_file_path = 'Neck/Neck_cross_section_u_at_0_percentage_worldspace.dat'
save_file_path = dirPath + 'Chest/Chest_cross_section_u_at_100_percentage_worldspace.dat'
source_file_path = 'Belly/Belly_cross_section_u_at_100_percentage_worldspace.dat'
cut_points = [18, 28, 52, 62]
angles = []
center = om.MVector()
Ca = []  # a parameters for x-axis, for every segment in a certain cross-section curve
Cb = []  # b parameters for y-axis, for every segment in a certain cross-section curve
Cc = []  # c parameters for z-axiz, for every segment in a certain cross-section curve

with open(dirPath+source_file_path, 'r') as f:
    line = f.readline()
    if not line.startswith('range:'):
        raise TypeError("bad format new segment should begin with range:")
    line = f.readline()
    if line.startswith('center:'):
        list = line.split(' ')
        center_x = float(list[-3])
        center_y = float(list[-2])
        center_z = float(list[-1])
        center = om.MVector(center_x, center_y, center_z)
        line = f.readline()

    if line.startswith('a:'):
        list = line.split(' ')
        for a in list[1:-1]:
            Ca.append(float(a))
        line = f.readline()
    else:
        raise ValueError("Bad file format, missing /'a:/' clause after the /'range:/' clause !")
    if line.startswith('b:'):
        list = line.split(' ')
        for b in list[1:-1]:
            Cb.append(float(b))
        line = f.readline()
    else:
        raise ValueError("Bad file format, missing /'b:/' clause after the /'a:/' clause !")
    if line.startswith('c:'):
        list = line.split(' ')
        for c in list[1:-1]:
            Cc.append(float(c))
        line = f.readline()
    else:
        raise ValueError("Bad file format, missing /'c:/' clause after the /'b:/' clause !")
    if line.startswith('angles:'):
        list = line.split(' ')
        for angle in list[1:-1]:
            angles.append(float(angle))
    else:
        raise ValueError("Bad file format, missing /'angles:/' clause after the /'c:/' clause !")

vertices_matrix, angles_matrix, segment_center_list = curve_fitting.cut_curve(vertices, angles, cut_points, isClosed = True)


#curve_fitting.maya_polygon_plane(composite_vertices_n)
for i in range(len(cut_points)):
    if i == 0:
        f = open(save_file_path, "w+")
    else:
        f = open(save_file_path, "a+")

    f.write('range:')
    f.write('{} to {} \n'.format(angles_matrix[i][0],angles_matrix[i][-1]))
    f.write('center: {} {} {}\n'.format(segment_center_list[i][0], segment_center_list[i][1], segment_center_list[i][2]))
    f.write('a: ')
    f.write(str(center[0] + Ca[0] - segment_center_list[i][0]) + ' ')
    for j in range(1, len(Ca)):
        f.write(str(Ca[j]) + ' ')
    f.write('\n')
    f.write('b: ')
    f.write(str(center[1] + Cb[0] - segment_center_list[i][1]) + ' ')
    for j in range(1, len(Cb)):
        f.write(str(Cb[j]) + ' ')
    f.write('\n')
    f.write('c: ')
    f.write(str(center[2] + Cc[0] - segment_center_list[i][2]) + ' ')
    for j in range(1, len(Cc)):
        f.write(str(Cc[j]) + ' ')
    f.write('\n')
    f.write('angles: ')
    for angle in angles_matrix[i]:
        f.write(str(angle) + ' ')
    f.write('\n')

    f.close()
