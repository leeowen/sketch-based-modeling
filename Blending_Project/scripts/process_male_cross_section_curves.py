import curve_fitting
#import maya.api.OpenMaya as om


def maya_polygon_plane(generalisedEllipseVertices):
    """
    test the plane in Maya
    """
    width = 10.0
    length = 10.0
    face_count = 1
    vertex_count = len(generalisedEllipseVertices)
    # Create vertex positions
    vertices = om.MFloatPointArray()
    for v in generalisedEllipseVertices:
        vertices.append(om.MFloatPoint(v[0] * 300, v[1] * 300, v[2] * 300))

    # Vertex count for this polygon face
    face_vertexes = om.MIntArray()
    face_vertexes.append(vertex_count)

    # Vertex indexes for this polygon face
    vertex_indexes = om.MIntArray()
    vertex_indexes.copy([i for i in range(vertex_count)])

    # Create mesh
    mesh_object = om.MObject()
    mesh = om.MFnMesh()
    mesh_object = mesh.create(vertices, face_vertexes, vertex_indexes)
    mesh.updateSurface()
    # Assign default shading
    cmds.sets(mesh.name(), edit=True, forceElement="initialShadingGroup")


"""    
file_paths = [
'Source_Chest_cross_section_u_at_0_percentage_worldspace.dat',
'Source_Chest_cross_section_u_at_10_percentage_worldspace.dat',
'Source_Chest_cross_section_u_at_22_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightArm_cross_section_u_at_35_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_55_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_95_percentage_worldspace.dat',
    'Source_RightArm_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightForeArm_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_35_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_RightForeArm_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftLeg_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftLeg_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightThigh_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_RightThigh_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftThigh_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_LeftThigh_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftKnee_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftKnee_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightKnee_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightKnee_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightAnkle_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightAnkle_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftAnkle_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftAnkle_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftFoot_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftFoot_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightFoot_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightFoot_cross_section_u_at_90_percentage_worldspace.dat'
]   

file_paths = [
    'Source_Neck_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_Neck_cross_section_u_at_98_percentage_worldspace.dat'
]

file_paths = [
    'Source_Hip_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_25_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_75_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_Belly_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_18_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_25_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftForeArm_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_35_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_LeftForeArm_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_LeftArm_cross_section_u_at_35_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_55_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_LeftArm_cross_section_u_at_100_percentage_worldspace.dat'
]

file_paths = [
    'Source_RightLeg_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_RightLeg_cross_section_u_at_90_percentage_worldspace.dat'
]
"""

file_paths = [
'Source_Chest_cross_section_u_at_22_percentage_worldspace.dat',
'Source_Chest_cross_section_u_at_66_percentage_worldspace.dat',
'Source_Chest_cross_section_u_at_86_percentage_worldspace.dat',
'Source_Chest_cross_section_u_at_96_percentage_worldspace.dat'
]


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

    file_name = file_path.replace('Source_', '')
    directory = dirPath + file_name.split('_cross_section_')[0]

    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    save_file_path = directory + '/' + file_name

    J = 6

    if 'worldspace' not in file_path:
        center = curve_fitting.getCenter(vertices)
        d_bar = curve_fitting.get_d_bar(vertices, center)
        angles = curve_fitting.calculateAngle_2D(vertices, center)
        a, b = curve_fitting.getCoefficients_2D(J, vertices, center, angles)
        generalisedEllipseVertices, Ea, Em = curve_fitting.formGeneralizedEllipse_2D(a, b, vertices, center, angles, d_bar, 0)

        with open(save_file_path, "w+") as f:
            f.write('range:')
            f.write('0-360 \n')
            f.write('center: {} {}\n'.format(center[0], center[2]))
            f.write('a: ')
            for i in a:
                f.write(str(i) + ' ')
            f.write('\n')
            f.write('b: ')
            for i in b:
                f.write(str(i) + ' ')
            f.write('\n')
            f.write('angles: ')
            for i in angles:
                f.write(str(i)+' ')
            f.write('\n')
    else:
        if 'Chest' in file_path:
            # delete the unnecessary points
            delete_points_list = {'Source_Chest_cross_section_u_at_22_percentage_worldspace.dat': [34, 37, 83, 86],
                                  'Source_Chest_cross_section_u_at_66_percentage_worldspace.dat': [30, 39, 81, 90],
                                  'Source_Chest_cross_section_u_at_86_percentage_worldspace.dat': [28, 40, 80, 92],
                                  'Source_Chest_cross_section_u_at_96_percentage_worldspace.dat': [28, 40, 80, 92]
                                  }
            Ea_criteria = 0.01
            Em_criteria = 0.025
            dir_path = cmds.workspace(fn=True) + '/data/'
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
                del vertices[delete_points_list.get(file_path)[2] + 1:delete_points_list.get(file_path)[3]]
                del vertices[delete_points_list.get(file_path)[0] + 1:delete_points_list.get(file_path)[1]]

                with open(dir_path + 'removed_' + file_path, "w") as f:
                    for v in vertices:
                        f.write('{} {} {}\n'.format(v[0], v[1], v[2]))
                    # rebuild the curve for every cross-section
                    center = curve_fitting.getCenter_3D(vertices)
                    angles = curve_fitting.calculateAngle_3D(vertices, center)
                    d_bar = curve_fitting.get_d_bar_3D(vertices, center)
                    J = [20, 2, 20]
                    coe = curve_fitting.getCoefficients_3D(J, vertices, center, angles)
                    # add points to replace the deleted points
                    delta_angle = (angles[delete_points_list.get(file_path)[0] + 1] - angles[delete_points_list.get(file_path)[0]]) / (delete_points_list.get(file_path)[1]-delete_points_list.get(file_path)[0])
                    for index in range(delete_points_list.get(file_path)[0],
                                       delete_points_list.get(file_path)[1] - 1):
                        angles.insert(index + 1, delta_angle + angles[index])
                    for index in range(delete_points_list.get(file_path)[2],
                                       delete_points_list.get(file_path)[3] - 1):
                        angles.insert(index + 1, delta_angle + angles[index])

                    tmp_x = curve_fitting.form_vertices_of_fragment_single(coe[0], center, angles, 0)
                    tmp_y = curve_fitting.form_vertices_of_fragment_single(coe[1], center, angles, 1)
                    tmp_z = curve_fitting.form_vertices_of_fragment_single(coe[2], center, angles, 2)
                    new_vertices = []
                    for i in range(len(tmp_x)):
                        new_vertices.append(om.MVector(tmp_x[i], tmp_y[i], tmp_z[i]))

                isClosed = True
                cut_points = [28, 40, 80, 92]
                # split data for N segments
                N = len(cut_points)
                vertices_matrix = []
                angles_matrix = []
                segment_center_list = []
                Ea_criteria = 0.01
                Em_criteria = 0.025
                # split curve
                for i in range(N - 1):
                    tmp_vertices = new_vertices[cut_points[i]:cut_points[(i + 1) % N] + 1]
                    tmp_angles = angles[cut_points[i]:cut_points[(i + 1) % N] + 1]
    
                    vertices_matrix.append(tmp_vertices)
                    angles_matrix.append(tmp_angles)
                    tmp_center = curve_fitting.getCenter_3D(tmp_vertices)
                    segment_center_list.append(tmp_center)
    
                tmp_vertices = new_vertices[cut_points[N - 1]:numPt]
                tmp_vertices.extend(new_vertices[0:cut_points[0] + 1])
                tmp_angles = angles[cut_points[N - 1]:numPt]
                tmp_angles.extend(angles[0:cut_points[0] + 1])
                vertices_matrix.append(tmp_vertices)
                angles_matrix.append(tmp_angles)
                tmp_center = curve_fitting.getCenter_3D(tmp_vertices)
                segment_center_list.append(tmp_center)
                d_bar = curve_fitting.get_d_bar_3D(new_vertices, center)

   
                # for the first segment
                composite_vertices = []
                coefficients = []
                getCoefficients_for_first_generalized_elliptic_segment = curve_fitting.getCoefficients_single
                J1, coe = curve_fitting.findJ_3D(vertices_matrix[0], angles_matrix[0], d_bar, segment_center_list[0], Ea_criteria, Em_criteria, cut_points[0])
                composite_vertices_0, Ea0, Em0 = curve_fitting.form_vertices_of_fragment_3D(coe, vertices_matrix[0], segment_center_list[0], angles_matrix[0], d_bar, cut_points[0])
                composite_vertices.append(composite_vertices_0)

                coefficients.append(coe)
                Ea = Ea + Ea0
                Em = Em + Em0
                print Ea, Em
                maya_polygon_plane(composite_vertices_0)
"""
                # for the segment(s) in-between
                for i in range(1, len(cut_points) - 1):
                    [x0_tan, y0_tan, z0_tan], [x0, y0, z0] = curve_fitting.position_and_tangent_of_parametric_point_3D(coefficients[i - 1], angles_matrix[i - 1][-1])
                    x0 += segment_center_list[i - 1][0]
                    y0 += segment_center_list[i - 1][1]
                    z0 += segment_center_list[i - 1][2]
                    cut_pt_index = cut_points[i - 1]
                    if isClosed == True:
                        cut_pt_index = cut_points[i]
                    previous = {'position x': x0, 'position y': y0, 'position z': z0, 'tangent x': x0_tan,
                                'tangent y': y0_tan, 'tangent z': z0_tan, 'cut point index': cut_pt_index}
                    J = curve_fitting.findJ_for_non_end_composite_3D(vertices_matrix[i], angles_matrix[i], d_bar, segment_center_list[i],
                                                     Ea_criteria, Em_criteria, previous, curve_fitting.getCoefficients_for_non_end_composite_single)
                    coe = curve_fitting.getCoefficients_for_non_end_composite_3D(J, vertices_matrix[i],
                                                                     segment_center_list[i], angles_matrix[i], previous)
                    composite_vertices_n, Ean, Emn = curve_fitting.form_vertices_of_fragment_3D(coe, vertices_matrix[i],
                                                                 segment_center_list[i], angles_matrix[i],d_bar, cut_pt_index)
                    coefficients.append(coe)
                    composite_vertices.append(composite_vertices_n)

                    previous_count = 0
                    for j in range(i):
                        previous_count += len(vertices_matrix[j])
                    Ea = (Ea * previous_count + Ean * len(vertices_matrix[-1])) / numPt
                    if Em < Emn:
                        Em = Emn
    
                    #maya_polygon_plane(composite_vertices_n)
                    #maya_polygon_plane(vertices_matrix[i])
    
                # for the end segment that links the first segment
                tan0, pos0 = position_and_tangent_of_parametric_point_3D(coefficients, angles_matrix[-2][-1])
                tan0[0] += segment_center_list[-2][0]
                tan0[1] += segment_center_list[-2][1]
                tan0[2] += segment_center_list[-2][2]
                previous = {'position x': pos0[0], 'position y': pos0[1], 'position z': pos0[2], 'tangent x': tan0[0], 'tangent y': tan0[1],
                            'tangent z': tan0[2], 'cut point index': cut_points[-1]}
                tan1, pos1 = position_and_tangent_of_parametric_point_3D(coefficients, angles_matrix[0][0])
                tan1[0] += segment_center_list[0][0]
                tan1[1] += segment_center_list[0][1]
                tan1[2] += segment_center_list[0][2]
                next = {'position x': pos1[0], 'position y': pos1[1], 'position z': pos1[2], 'tangent x': tan1[0], 'tangent y': tan1[1],
                        'tangent z': tan1[2], 'cut point index': cut_points[0]}
    
                Jn = findJ_for_end_segment_3D(vertices_matrix[-1], angles_matrix[-1], d_bar, segment_center_list[-1], Ea_criteria, Em_criteria, previous, next)
    
                coe = getCoefficients_for_end_composite(Jn, vertices_matrix[-1], segment_center_list[-1],angles_matrix[-1], previous, next)
                composite_vertices_n, Ean, Emn = form_vertices_of_fragment_3D(coe, vertices_matrix[-1], segment_center_list[-1], angles_matrix[-1], d_bar,cut_points[-1])
                coefficients.append(coe)
                composite_vertices.append(composite_vertices_n)
    
                Ea = (Ea * (numPt - len(vertices_matrix[-1])) + Ean * len(vertices_matrix[-1])) / numPt
                if Em < Emn:
                    Em = Emn
                manualJ_value = J_total + Jn

                with open(save_file_path, "w+") as f:
                    for i in range(N - 1):
                        f.write('range:')
                        f.write('{} to {} \n'.format(angles_matrix[i][0], angles_matrix[i][-1]))
                        f.write('center: {} {} {}\n'.format(segment_center_list[i][0], segment_center_list[i][1], segment_center_list[i][2]))
                        f.write('a: ')
                        for j in coefficients[i][0]:
                            f.write(str(j) + ' ')
                        f.write('\n')
                        f.write('b: ')
                        for j in coefficients[i][1]:
                            f.write(str(j) + ' ')
                        f.write('\n')
                        f.write('c: ')
                        for j in coefficients[i][2]:
                            f.write(str(j) + ' ')
                        f.write('\n')
                        f.write('angles: ')
                        for j in angles_matrix[i]:
                           f.write(str(j)+' ')
                        f.write('\n')
        else:
            if 'Foot' in file_path:
                J = 2
                coe = curve_fitting.getCoefficients3_swap_yz(J, vertices, center, angles)
                generalisedEllipseVertices, Ea, Em = curve_fitting.formGeneralizedEllipse3_swap_yz(coe, vertices, center, angles, d_bar)
            else:
                if 'Neck' in file_path:
                    J = [8,2,8]
                if 'Hip' or 'Belly' in file_path:
                    J = [18,2,18]
                if 'Arm' in file_path:
                    J = [8,2,8]
                coe = curve_fitting.getCoefficients_3D(J, vertices, center, angles)
                generalisedEllipseVertices, Ea, Em = curve_fitting.formGeneralizedEllipse_3D(coe, vertices, center, angles, d_bar)

            with open(save_file_path, "w+") as f:
                f.write('range:')
                f.write('0-360 \n')
                f.write('center: {} {} {}\n'.format(center[0], center[1], center[2]))
                f.write('a: ')
                for i in coe[0]:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('b: ')
                for i in coe[1]:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('c: ')
                for i in coe[2]:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('angles: ')
                for i in angles:
                   f.write(str(i)+' ')
                f.write('\n')
"""

