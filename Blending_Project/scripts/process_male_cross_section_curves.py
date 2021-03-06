import maya.api.OpenMaya as om
import os
sys.path.append(cmds.workspace(fn=True)+'/scripts/')
import curve_fitting
import math


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
    'Source_Neck_cross_section_u_at_100_percentage_worldspace.dat'
]


file_paths = [
    'Source_Shoulder_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_Shoulder_cross_section_u_at_90_percentage_worldspace.dat'
]

file_paths = [
    'Source_Hip_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_70_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_80_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_90_percentage_worldspace.dat',
    'Source_Hip_cross_section_u_at_100_percentage_worldspace.dat'
]
"""
file_paths = [
    'Source_Belly_cross_section_u_at_0_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_10_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_20_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_30_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_40_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_50_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_60_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_68_percentage_worldspace.dat',
    'Source_Belly_cross_section_u_at_75_percentage_worldspace.dat'
]
"""
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

file_paths = [
'modified_Source_Chest_cross_section_u_at_0_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_22_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_66_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_86_percentage_worldspace.dat',
'modified_Source_Chest_cross_section_u_at_100_percentage_worldspace.dat'
]
"""
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
    directory = dirPath + file_name.split('_cross_section_')[0]

    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    save_file_path = directory + '/' + file_name

    if 'worldspace' not in file_path:
        J = 6
        center = curve_fitting.getCenter(vertices)
        d_bar = curve_fitting.get_d_bar(vertices, center)
        #angles = curve_fitting.calculateAngle_2D(vertices, center)
        angles = [i*2*math.pi/len(vertices) for i in range(len(vertices))]
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
            cut_points = [27, 42, 78, 93]
            isClosed = True

            #curve_fitting.maya_polygon_plane(vertices)

            # cut curve
            center = curve_fitting.getCenter_3D(vertices)
            angles = [i*2*math.pi/numPt for i in range(numPt)]
            #angles = curve_fitting.calculateAngle_3D(vertices, center)
            d_bar = curve_fitting.get_d_bar_3D(vertices, center)
            vertices_matrix, angles_matrix, segment_center_list = curve_fitting.cut_curve(vertices, angles, cut_points, isClosed)

            composite_vertices = []
            composite_coefficients = []


            # for the first segment
            J0, coe = curve_fitting.findJ_3D(vertices_matrix[0], angles_matrix[0], d_bar, segment_center_list[0], Ea_criteria, Em_criteria, cut_points[0])
            composite_vertices_0 = curve_fitting.form_vertices_of_fragment_3D(coe, segment_center_list[0], angles_matrix[0])
            Ea0, Em0 = curve_fitting.calculate_Ea_Em_3D(vertices_matrix[0], composite_vertices_0, d_bar, cut_points[0])
            composite_vertices.append(composite_vertices_0)
            composite_coefficients.append(coe)
            Ea = Ea + Ea0
            Em = Em + Em0
            curve_fitting.maya_polygon_plane(composite_vertices_0)
            print "J0 = {}".format(J0)

            # for in-between segment(s)
            for i in range(1, len(cut_points) - 1):
                cut_pt_index = cut_points[i - 1]
                tan0, p0 = curve_fitting.position_and_tangent_of_parametric_point_3D(composite_coefficients[-1],
                                                                       angles_matrix[i - 1][-1])
                p0[0] += segment_center_list[i - 1][0]
                p0[1] += segment_center_list[i - 1][1]
                p0[2] += segment_center_list[i - 1][2]

                cut_pt_index = cut_points[i]
                previous = {'position x': p0[0], 'position y': p0[1], 'position z': p0[2], 'tangent x': tan0[0],
                            'tangent y': tan0[1], 'tangent z': tan0[2],
                            'cut point index': cut_pt_index}
                J, coen = curve_fitting.findJ_for_non_end_composite_3D(vertices_matrix[i], angles_matrix[i], d_bar, segment_center_list[i], Ea_criteria, Em_criteria, previous, curve_fitting.getCoefficients_for_non_end_composite_single)
                composite_vertices_n = curve_fitting.form_vertices_of_fragment_3D(coen, segment_center_list[i], angles_matrix[i])
                Ean, Emn = curve_fitting.calculate_Ea_Em_3D(vertices_matrix[i], composite_vertices_n, d_bar, cut_points[i])
                Ea = (Ea * cut_pt_index + Ean * len(composite_vertices_n)) / (cut_pt_index + len(composite_vertices_n))
                Em = max(Em, Emn)
                composite_vertices.append(composite_vertices_n)
                composite_coefficients.append(coen)
                print "J{} = {}".format(i, J)
                curve_fitting.maya_polygon_plane(composite_vertices_n)

            # for the end segment that links the first segment
            tan0, p0 = curve_fitting.position_and_tangent_of_parametric_point_3D(composite_coefficients[-1], angles_matrix[-2][-1])
            p0[0] += segment_center_list[-2][0]
            p0[1] += segment_center_list[-2][1]
            p0[2] += segment_center_list[-2][2]
            previous = {'position x': p0[0], 'position y': p0[1], 'position z': p0[2], 'tangent x': tan0[0], 'tangent y': tan0[1], 'tangent z': tan0[2], 'cut point index': cut_points[-1]}
            tan1, p1 = curve_fitting.position_and_tangent_of_parametric_point_3D(composite_coefficients[0], angles_matrix[0][0])
            p1[0] += segment_center_list[0][0]
            p1[1] += segment_center_list[0][1]
            p1[2] += segment_center_list[0][2]
            next = {'position x': p1[0], 'position y': p1[1], 'position z': p1[2], 'tangent x': tan1[0], 'tangent y': tan1[1], 'tangent z': tan1[2], 'cut point index': cut_points[0]}

            Jn = curve_fitting.findJ_for_end_segment_3D(vertices_matrix[-1], angles_matrix[-1], d_bar, segment_center_list[-1], Ea_criteria, Em_criteria, previous, next)
            print "J{} = {}".format(len(cut_points) - 1, Jn)
            coen = curve_fitting.getCoefficients_for_end_composite_3D(Jn, vertices_matrix[-1], segment_center_list[-1], angles_matrix[-1], previous, next)
            composite_vertices_n = curve_fitting.form_vertices_of_fragment_3D(coen, segment_center_list[-1], angles_matrix[-1])
            Ean, Emn = curve_fitting. calculate_Ea_Em_3D(vertices_matrix[-1], composite_vertices_n, d_bar, cut_points[-1])
            composite_coefficients.append(coen)
            composite_vertices.append(composite_vertices_n)

            Ea = (Ea * (numPt - len(vertices_matrix[-1])) + Ean * len(vertices_matrix[-1])) / numPt
            if Em < Emn:
                Em = Emn

            curve_fitting.maya_polygon_plane(composite_vertices_n)
            for i in range(len(cut_points)):
                if i == 0:
                    f = open(save_file_path, "w+")
                else:
                    f = open(save_file_path, "a+")

                f.write('range:')
                f.write('{} to {} \n'.format(angles_matrix[i][0],angles_matrix[i][-1]))
                f.write('center: {} {} {}\n'.format(segment_center_list[i][0], segment_center_list[i][1], segment_center_list[i][2]))
                f.write('a: ')
                for a in composite_coefficients[i][0]:
                    f.write(str(a) + ' ')
                f.write('\n')
                f.write('b: ')
                for b in composite_coefficients[i][1]:
                    f.write(str(b) + ' ')
                f.write('\n')
                f.write('c: ')
                for c in composite_coefficients[i][2]:
                    f.write(str(c) + ' ')
                f.write('\n')
                f.write('angles: ')
                for angle in angles_matrix[i]:
                    f.write(str(angle) + ' ')
                f.write('\n')

                f.close()
        else:
            isClosed = True
            # curve_fitting.maya_polygon_plane(vertices)

            # cut curve
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
