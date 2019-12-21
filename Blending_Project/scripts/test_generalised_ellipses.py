import curve_fitting
import maya.api.OpenMaya as OpenMaya

def maya_polygon_plane(generalisedEllipseVertices):
    width = 10.0
    length = 10.0
    face_count = 1
    vertex_count = len(generalisedEllipseVertices)
    # Create vertex positions
    vertices = OpenMaya.MFloatPointArray()
    for v in generalisedEllipseVertices:
        vertices.append(OpenMaya.MFloatPoint(v[0]*300, 0.0 , v[1]*300))

    # Vertex count for this polygon face
    face_vertexes = OpenMaya.MIntArray()
    face_vertexes.append(vertex_count)

    # Vertex indexes for this polygon face
    vertex_indexes = OpenMaya.MIntArray()
    vertex_indexes.copy([i for i in range(vertex_count)])

    # Create mesh
    mesh_object = OpenMaya.MObject()
    mesh = OpenMaya.MFnMesh()
    mesh_object = mesh.create(vertices, face_vertexes, vertex_indexes)
    mesh.updateSurface()
    # Assign default shading
    cmds.sets(mesh.name(), edit=True, forceElement="initialShadingGroup")


if __name__ == "__main__":
    """
    file_paths = [
        'Source_LeftForeArm_cross_section_u_at_0_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_10_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_20_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_30_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_35_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_40_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_50_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_60_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_70_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_80_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_90_percentage.dat',
        'Source_LeftForeArm_cross_section_u_at_100_percentage.dat'
    ]
    """
    file_paths = [
        'Source_LeftLeg_cross_section_u_at_10_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_20_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_30_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_40_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_50_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_60_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_70_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_80_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_90_percentage.dat',
        'Source_LeftLeg_cross_section_u_at_100_percentage.dat'
    ]
    """
    file_paths = [
        'Source_LeftArm_cross_section_u_at_16_percentage.dat',
        'Source_LeftArm_cross_section_u_at_20_percentage.dat',
        'Source_LeftArm_cross_section_u_at_30_percentage.dat',
        'Source_LeftArm_cross_section_u_at_40_percentage.dat',
        'Source_LeftArm_cross_section_u_at_51_percentage.dat',
        'Source_LeftArm_cross_section_u_at_60_percentage.dat',
        'Source_LeftArm_cross_section_u_at_70_percentage.dat',
        'Source_LeftArm_cross_section_u_at_80_percentage.dat',
        'Source_LeftArm_cross_section_u_at_90_percentage.dat',
        'Source_LeftArm_cross_section_u_at_100_percentage.dat'
    ]
    """
    dirPath = cmds.workspace(fn=True)+'/data/'
    for file_path in file_paths:
        vertices = []
        numPt = 0
        with open(dirPath + file_path, 'r') as f:
            content = f.readlines()
            for line in content:
                p = line.split()
                numPt += 1
                vertices.append(om.MVector(float(p[0]), float(p[1]), float(p[2])))
    
            center = curve_fitting.getCenter(vertices)
            d_bar = curve_fitting.get_d_bar(vertices, center)
            angles = curve_fitting.calculateAngle(vertices, center)
            J = 6
            a, b = curve_fitting.getCoefficients(J, vertices, center, angles)
            generalisedEllipseVertices, Ea, Em = curve_fitting.formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)

            maya_polygon_plane(generalisedEllipseVertices)