import math,sys, os
import operator
from PySide2 import QtCore, QtGui
import maya.api.OpenMaya as om
import maya.cmds as cmds

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def maya_useNewAPI():
    pass


class IllegalArgumentError(ValueError):
    pass


def getCenter_2D(vertices):
    center = om.MVector(0.0, 0.0, 0.0)
    for v in vertices:
        center += om.MVector(float(v[0]), 0, float(v[2]))
    center = center / len(vertices)
    return center


def getCenter_3D(vertices):
    center = om.MVector(0.0, 0.0, 0.0)
    for v in vertices:
        center = center + v
    center = center / len(vertices)
    return center


def get_d_bar_2D(vertices,center):
    d_bar = []
    for v in vertices:
        d_bar.append(math.sqrt((v[0] - center[0]) ** 2 + (v[2] - center[2]) ** 2))
    return d_bar


def get_d_bar_3D(vertices,center):
    d_bar = []
    for v in vertices:
        d_bar.append(math.sqrt((v[0] - center[0]) ** 2 + (v[1] - center[1]) ** 2 + (v[2] - center[2]) ** 2))
    return d_bar


def calculateAngle_2D(vertices,center):
    angles = []
    for v in vertices:
        anglem = (v[2] - center[2]) / math.sqrt(
            (v[2] - center[2]) ** 2 + (v[0] - center[0]) ** 2)
        anglem = math.acos(anglem)

        if v[0] > center[0] and v[2] > center[2]:
            pass
        elif v[0] > center[0] and v[2] < center[2]:
            pass
        elif v[0] < center[0] and v[2] > center[2]:
            anglem = 2 * math.pi - anglem
        elif v[0] < center[0] and v[2] < center[2]:
            anglem = 2 * math.pi - anglem

        while anglem > 2 * math.pi:
            anglem = anglem - 2 * math.pi
        while anglem < 0:
            anglem = anglem + 2 * math.pi

        angles.append(anglem)
    return angles


def calculateAngle_3D(vertices,center):
    return calculateAngle_2D(vertices,center)


# need to be improved
def calculateArcLength_3D(vertices):
    totalArcLength = 0.0
    arcLength = []
    numPt = len(vertices)
    for i in range(0, numPt):
        arcL = math.sqrt(pow(vertices[(i + 1) % numPt][0] - vertices[i][0], 2) + pow(
            vertices[(i + 1) % numPt][2] - vertices[i][2], 2))
        arcLength.append(arcL)
        totalArcLength += arcLength
    return totalArcLength, arcLength


def calculateTangent_2D(vertices, angles):
    tangents = []
    numPt = len(vertices)
    for i in range(0, numPt):
        # primes refer to the derivatives dx/dt, dy/dt with respect to the angle t
        dt = angles[(i + 1) % numPt] - angles[i - 1]
        if dt < 0:
            dt += math.pi*2
        tan_x = (vertices[(i + 1) % numPt][0] - vertices[i - 1][0]) / dt
        tan_y = (vertices[(i + 1) % numPt][2] - vertices[i - 1][2]) / dt
        tangents.append((tan_x, tan_y))
    return tangents


def calculateNormal_2D(tangents):
    normals=[]
    for t in tangents:
        normal = QtGui.QVector2D(-t[1], t[0]).normalized()
        normals.append(normal)
    return normals


def calculateCurvature_2D(tangents, angles):
    curvatures = []
    numPt = len(tangents)
    for i in range(0, numPt):
        dt = angles[(i + 1) % numPt] - angles[i - 1]
        if dt < 0:
            dt += math.pi*2
        d2xdt2 = (tangents[(i + 1) % numPt][0] - tangents[i][0]) / dt
        d2ydt2 = (tangents[(i + 1) % numPt][1] - tangents[i][1]) / dt
        k = (tangents[i][0] * d2ydt2 - tangents[i][1] * d2xdt2) / pow(
            pow(tangents[i][0], 2) + pow(tangents[i][1], 2), 3.0 / 2.0)
        curvatures.append(k)
    return curvatures


def sortCurvature(curvatures):
    dict = {}

    for i in range(len(curvatures)):
        dict[i] = abs(curvatures[i])

    sorted_d = sorted(dict.items(), key = operator.itemgetter(1))
    return sorted_d


def cut_curve(vertices, angles, cut_points, isClosed):
    N = len(cut_points)
    numPt = len(vertices)
    if isClosed == True:
        if N < 2:
            raise ValueError('for closed curve, at least 2 cut points are required')
    else:
        if N < 1:
            raise ValueError('for open curve, at least 1 cut point is required')
    cut_points.sort()  # ascend
    # split data for N segments

    vertices_matrix = []
    angles_matrix = []
    segment_center_list = []

    if isClosed == True:
        for i in range(N - 1):
            tmp_vertices = vertices[cut_points[i]:cut_points[(i + 1) % N] + 1]
            tmp_angles = angles[cut_points[i]:cut_points[(i + 1) % N] + 1]

            vertices_matrix.append(tmp_vertices)
            angles_matrix.append(tmp_angles)
            tmp_center = getCenter_3D(tmp_vertices)
            segment_center_list.append(tmp_center)

        tmp_vertices = vertices[cut_points[N - 1]:numPt]
        tmp_vertices.extend(vertices[0:cut_points[0] + 1])
        tmp_angles = angles[cut_points[N - 1]:numPt]
        tmp_angles.extend(angles[0:cut_points[0] + 1])
        vertices_matrix.append(tmp_vertices)
        angles_matrix.append(tmp_angles)
        tmp_center = getCenter_3D(tmp_vertices)
        segment_center_list.append(tmp_center)

    else:
        tmp_vertices = vertices[0:cut_points[0] + 1]
        tmp_angles = angles[0:cut_points[0] + 1]
        vertices_matrix.append(tmp_vertices)
        angles_matrix.append(tmp_angles)
        tmp_center = getCenter_3D(tmp_vertices)
        segment_center_list.append(tmp_center)

        if N > 1:
            for i in range(N - 1):
                tmp_vertices = vertices[cut_points[i]:cut_points[i + 1] + 1]
                tmp_angles = angles[cut_points[i]:cut_points[i + 1] + 1]
                vertices_matrix.append(tmp_vertices)
                angles_matrix.append(tmp_angles)
                tmp_center = getCenter_3D(tmp_vertices)
                segment_center_list.append(tmp_center)

        tmp_vertices = vertices[cut_points[N - 1]: numPt]
        tmp_angles = angles[cut_points[N - 1]: numPt]
        vertices_matrix.append(tmp_vertices)
        angles_matrix.append(tmp_angles)
        tmp_center = getCenter_3D(tmp_vertices)
        segment_center_list.append(tmp_center)

    return vertices_matrix, angles_matrix, segment_center_list


def formGeneralizedEllipse_2D(a, b, vertices, center, angles, d_bar,index=0):
    # CoefficientMatrix
    I = len(vertices)
    generalisedEllipseVertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    J = len(a) / 2
    for i in range(I):
        generalisedEllipseVertices[i][0] = center[0] + a[0]
        generalisedEllipseVertices[i][2] = center[2] + b[0]
        v = angles[i]
        for j in range(1, J + 1):
            generalisedEllipseVertices[i][0] += a[2 * j - 1] * math.cos(j * v) + a[2 * j] * math.sin(j * v)
            generalisedEllipseVertices[i][2] += b[2 * j - 1] * math.sin(j * v) + b[2 * j] * math.cos(j * v)

        di = math.sqrt((vertices[i][0] - generalisedEllipseVertices[i][0]) ** 2 + (
                    vertices[i][2] - generalisedEllipseVertices[i][2]) ** 2)
        Ea += (di / d_bar[(i+index) % I])
        if Em < di / d_bar[(i+index) % I]:
            Em = di / d_bar[(i+index) % I]
    Ea = Ea / I
    return generalisedEllipseVertices, Ea, Em


def formGeneralizedEllipse_3D(coe, vertices, center, angles, d_bar):
    # CoefficientMatrix
    I = len(vertices)
    generalisedEllipseVertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    J = len(a) / 2
    Jy = len(b) / 2
    for i in range(I):
        generalisedEllipseVertices[i][0] = center[0] + coe[0][0]
        generalisedEllipseVertices[i][1] = center[1] + coe[1][0]
        generalisedEllipseVertices[i][2] = center[2] + coe[2][0]
        v = angles[i]
        for j in range(1, J + 1):
            generalisedEllipseVertices[i][0] += a[2 * j - 1] * math.cos(j * v) + a[2 * j] * math.sin(j * v)
            generalisedEllipseVertices[i][2] += c[2 * j - 1] * math.cos(j * v) + c[2 * j] * math.sin(j * v)
        for j in range(1, Jy + 1):
            generalisedEllipseVertices[i][1] += b[2 * j - 1] * math.cos(j * v) + b[2 * j] * math.sin(j * v)

        di = math.sqrt((vertices[i][0] - generalisedEllipseVertices[i][0]) ** 2 + (
                    vertices[i][1] - generalisedEllipseVertices[i][1]) ** 2 + (
                    vertices[i][2] - generalisedEllipseVertices[i][2]) ** 2)
        Ea += (di / d_bar[i])
        if Em < di / d_bar[i]:
            Em = di / d_bar[i]
    Ea = Ea / I
    return generalisedEllipseVertices, Ea, Em


def formGeneralizedEllipse3_swap_yz(a, b, c, vertices, center, angles, d_bar):
    # CoefficientMatrix
    I = len(vertices)
    generalisedEllipseVertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    J = len(a) / 2
    Jz = len(c) / 2
    for i in range(I):
        generalisedEllipseVertices[i][0] = center[0] + a[0]
        generalisedEllipseVertices[i][1] = center[1] + b[0]
        generalisedEllipseVertices[i][2] = center[2] + c[0]
        v = angles[i]
        for j in range(1, J + 1):
            generalisedEllipseVertices[i][0] += a[2 * j - 1] * math.cos(j * v) + a[2 * j] * math.sin(j * v)
            generalisedEllipseVertices[i][1] += b[2 * j - 1] * math.cos(j * v) + b[2 * j] * math.sin(j * v)
        for j in range(1, Jz + 1):
            generalisedEllipseVertices[i][2] += c[2 * j - 1] * math.cos(j * v) + c[2 * j] * math.sin(j * v)

        di = math.sqrt((vertices[i][0] - generalisedEllipseVertices[i][0]) ** 2 + (
                    vertices[i][1] - generalisedEllipseVertices[i][1]) ** 2 + (
                    vertices[i][2] - generalisedEllipseVertices[i][2]) ** 2)
        d.append(di)
        Ea += (di / d_bar[i])
        if Em < di / d_bar[i]:
            Em = di / d_bar[i]
    Ea = Ea / I
    return generalisedEllipseVertices, Ea, Em


def extract_fragment_data(vertices, angles, fragment_range):
    numPt = len(vertices)
    start_index = int(numPt * fragment_range)
    end_index = start_index + numPt / 2

    if end_index >= numPt:
        end_index -= numPt

    # extract data
    angles_fragment = []
    vertices_fragment = []
    if start_index < end_index:
        for i in range(start_index, end_index + 1):
            angles_fragment.append(angles[i])
            vertices_fragment.append(vertices[i])
    else:
        for i in range(start_index, numPt):
            angles_fragment.append(angles[i])
            vertices_fragment.append(vertices[i])
        for i in range(end_index + 1):
            angles_fragment.append(angles[i])
            vertices_fragment.append(vertices[i])

    center_fragment = om.MVector(0.0, 0.0, 0.0)
    for v in vertices_fragment:
        center_fragment[0]= center_fragment[0] + v[0]
        center_fragment[2]= center_fragment[2] + v[2]
    center_fragment[0] = center_fragment[0] / len(vertices_fragment)
    center_fragment[2] = center_fragment[2] / len(vertices_fragment)

    return angles_fragment, vertices_fragment, center_fragment, start_index, end_index


def getCoefficients_2D(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    if J<3:
        raise IllegalArgumentError('J must be bigger than 3, you input {0}'.format(J))
    I = len(vertices)

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

    for i in range(I):
        bCoefficientMatrix[0, i] = 1.

    for i in range(I):  # for aCoefficientMatrix's column
        vi = angles[i]
        for j in range(1, J + 1):  # for aCoefficientMatrix's row
            bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
            bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][2] - center[2]) * math.sin(vi * j)
            bConstArray[2 * j] += (vertices[i][2] - center[2]) * math.cos(vi * j)

    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    a = getCoefficients_single([J,J], vertices, center, angles, 0)
    return a, b


def getCoefficients_3D(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    coe=[]
    for axis in [0,1,2]:# x,y,z axis
        coe.append(getCoefficients_single(J,vertices,center,angles,axis))
    return coe


def getCoefficients3_swap_yz(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    I = len(vertices)
    aConstArray = np.zeros(2 * J + 1)
    aCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

    Jz = 1
    cConstArray = np.zeros(2 * Jz + 1)
    cCoefficientMatrix = np.ndarray(shape=(2 * Jz + 1, I), dtype=float, order='C')

    for i in range(I):
        aCoefficientMatrix[0, i] = 1.
        bCoefficientMatrix[0, i] = 1.
        cCoefficientMatrix[0, i] = 1.

    for i in range(I):  # for aCoefficientMatrix's column
        vi = angles[i]
        for j in range(1, J + 1):  # for aCoefficientMatrix's row
            aCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            aCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
            aConstArray[2 * j - 1] += (vertices[i][0] - center[0]) * math.cos(vi * j)
            aConstArray[2 * j] += (vertices[i][0] - center[0]) * math.sin(vi * j)

            bCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            bCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][1] - center[1]) * math.cos(vi * j)
            bConstArray[2 * j] += (vertices[i][1] - center[1]) * math.sin(vi * j)
        for j in range(1, Jz + 1):
            cCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            cCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            cConstArray[2 * j - 1] += (vertices[i][2] - center[2]) * math.cos(vi * j)
            cConstArray[2 * j] += (vertices[i][2] - center[2]) * math.sin(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)
    C = np.dot(cCoefficientMatrix, cCoefficientMatrix.transpose())
    c = np.linalg.solve(C, cConstArray)

    return a, b, c


def getCoefficients_single(J, vertices, center, angles, axis):
    JJ = 0
    if isinstance(J,list):
        JJ = J[axis]
    elif isinstance(J,int):
        JJ = J

    if JJ < 0:
        raise IllegalArgumentError('J must be no smaller than 0, you input {0}'.format(JJ))
    I = len(vertices)
    ConstArray = np.zeros(2 * JJ + 1)
    CoefficientMatrix = np.ndarray(shape=(2 * JJ + 1, I), dtype=float, order='C')  # row-major

    for i in range(I):
        CoefficientMatrix[0, i] = 1.

    for i in range(I):  # for CoefficientMatrix's column
        vi = angles[i]
        for j in range(1, JJ + 1):  # for aCoefficientMatrix's row
            CoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            CoefficientMatrix[2 * j, i] = math.sin(vi * j)

            # ConstAtrray[0] always equals to 0 by definition!
            ConstArray[2 * j - 1] += (vertices[i][axis] - center[axis]) * math.cos(vi * j)
            ConstArray[2 * j] += (vertices[i][axis] - center[axis]) * math.sin(vi * j)

    M = np.dot(CoefficientMatrix, CoefficientMatrix.transpose())
    coe = np.linalg.solve(M, ConstArray)

    return coe


def getCoefficients_for_second_half_of_symmetrical_ellipse(a_first_half, b_first_half):
    a = np.zeros(len(a_first_half))
    b = np.zeros(len(b_first_half))
    a[0]=-a_first_half[0]
    b[0]=b_first_half[0]
    for j in range(1,len(a_first_half)/2+1):
        a[2*j-1]=-a_first_half[2*j-1]
        a[2*j]=a_first_half[2*j]
        b[2*j-1]=-b_first_half[2*j-1]
        b[2*j]=b_first_half[2*j]

    return a, b


# connected at v=0 and v=math.pi
def getCoefficients_for_second_half_of_composite_segment(J, vertices_second_half, angles_second_half, center_first_half, center_second_half, a_first_half, b_first_half):
    """
    the case when the second segment starts at v=math.pi and ends at v=0
    :param J:
    :param vertices_second_half:
    :param angles_second_half:
    :param center_first_half:
    :param center_second_half:
    :param a_first_half:
    :param b_first_half:
    :return:
    """
    I = len(vertices_second_half)
    J1 = (len(a_first_half)-1) / 2
    aConstArray = np.zeros(I)
    aCoefficientMatrix = np.ndarray(shape=(2 * J - 3, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(I)
    bCoefficientMatrix = np.ndarray(shape=(2 * J - 3, I), dtype=float, order='C')
    for i in range(len(vertices_second_half)):
        vi = angles_second_half[i]
        Cai = vertices_second_half[i][0] - (center_first_half[0] + a_first_half[0]) * math.cos(2 * vi) - center_second_half[0] * (1 - math.cos(2 * vi))
        Cbi = vertices_second_half[i][2] - (center_first_half[2] + b_first_half[0]) * math.cos(2 * vi) - center_second_half[2] * (1 - math.cos(2 * vi))
        aCoefficientMatrix[0, i] = 1 - math.cos(2 * vi)
        bCoefficientMatrix[0, i] = 1 - math.cos(2 * vi)
        for j in range(1, J1 + 1):
            D1 = (1 - (-1)**j) * math.cos(vi) + (1 + (-1)**j) * math.cos(2 * vi)
            D2 = (1 - (-1)**j) * math.sin(vi) + (1 + (-1)**j) * math.sin(2 * vi) / 2.0
            Cai = Cai - (D1 * a_first_half[2 * j - 1] + j * D2 * a_first_half[2 * j]) / 2.0
            Cbi = Cbi - (j * D2 * b_first_half[2 * j - 1] + D1 * b_first_half[2 * j]) / 2.0
        aConstArray[i] = Cai
        bConstArray[i] = Cbi
        for j in range(3, J + 1):
            D1 = (1 - (-1)**j) * math.cos(vi) + (1 + (-1)**j) * math.cos(2 * vi)
            D2 = (1 - (-1)**j) * math.sin(vi) + (1 + (-1)**j) * math.sin(2 * vi) / 2.0
            aCoefficientMatrix[2 * j - 5, i] = math.cos(j * vi) - D1 / 2.0
            aCoefficientMatrix[2 * j - 4, i] = math.sin(j * vi) - D2 * j / 2.0
            bCoefficientMatrix[2 * j - 5, i] = math.sin(j * vi) - D2 * j / 2.0
            bCoefficientMatrix[2 * j - 4, i] = math.cos(j * vi) - D1 / 2.0

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, np.dot(aCoefficientMatrix, aConstArray))
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, np.dot(bCoefficientMatrix, bConstArray))

    a1, a2, a4, b1, b2, b3 = 0, 0, 0, 0, 0, 0
    a3 = center_first_half[0] - center_second_half[0] + a_first_half[0] - a[0]
    b4 = center_first_half[2] - center_second_half[2] + b_first_half[0] - b[0]
    for j in range(1, J1 + 1):
        a1 = a1 + (1 - (-1)**j) * a_first_half[2 * j - 1] / 2.0
        a2 = a2 + (1 - (-1)**j) * j * a_first_half[2 * j] / 2.0
        a3 = a3 + (1 + (-1)**j) * a_first_half[2 * j - 1] / 2.0
        a4 = a4 + (1 + (-1)**j) * j * a_first_half[2 * j] / 4.0
        b1 = b1 + (1 - (-1)**j) * j * b_first_half[2 * j - 1] / 2.0
        b2 = b2 + (1 - (-1)**j) * b_first_half[2 * j] / 2.0
        b3 = b3 + (1 + (-1)**j) * j * b_first_half[2 * j - 1] / 4.0
        b4 = b4 + (1 + (-1)**j) * b_first_half[2 * j] / 2.0
    for j in range(3, J + 1):
        a1 = a1 - (1 - (-1)**j) * a[2 * j - 5] / 2.0
        a2 = a2 - (1 - (-1)**j) * j * a[2 * j - 4] / 2.0
        a3 = a3 - (1 + (-1)**j) * a[2 * j - 5] / 2.0
        a4 = a4 - (1 + (-1)**j) * j * a[2 * j - 4] / 4.0
        b1 = b1 - (1 - (-1)**j) * j * b[2 * j - 5] / 2.0
        b2 = b2 - (1 - (-1)**j) * b[2 * j - 4] / 2.0
        b3 = b3 - (1 + (-1)**j) * j * b[2 * j - 5] / 4.0
        b4 = b4 - (1 + (-1)**j) * b[2 * j - 4] / 2.0

    a = np.insert(a, 1, a1)
    a = np.insert(a, 2, a2)
    a = np.insert(a, 3, a3)
    a = np.insert(a, 4, a4)
    b = np.insert(b, 1, b1)
    b = np.insert(b, 2, b2)
    b = np.insert(b, 3, b3)
    b = np.insert(b, 4, b4)

    return a, b


def getCoefficients_for_fragmented_ellipse( J, angles, vertices, center):
    I = len(angles)
    aConstArray = np.zeros(2 * J + 1)
    aCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

    for i in range(I):
        aCoefficientMatrix[0, i] = 1.
        bCoefficientMatrix[0, i] = 1.

    for i in range(I):  # for aCoefficientMatrix's column
        for j in range(1, J + 1):  # for aCoefficientMatrix's row
            vi = angles[i]

            aCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            aCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
            aConstArray[2 * j - 1] += (vertices[i][0] - center[0]) * math.cos(vi * j)
            aConstArray[2 * j] += (vertices[i][0] - center[0]) * math.sin(vi * j)

            bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
            bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][2] - center[2]) * math.sin(vi * j)
            bConstArray[2 * j] += (vertices[i][2] - center[2]) * math.cos(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    return a, b


def getCoefficients_for_end_composite_2D(J, vertices, center, angles, previous, next):
    I = len(vertices)
    e = 1e-10  # a value that is very close to zero

    def fun_y(y, *args):  # function to be minimized
        Ey = 0
        vertices = args[0]
        center = args[1]
        angles = args[2]
        J = args[3]
        for i in range(1, I - 1):
            yi = center[2] + y[0]
            for j in range(1, J + 1):
                yi += (y[2 * j - 1] * math.sin(j * angles[i]) + y[2 * j] * math.cos(j * angles[i]))
            Ey += (vertices[i][2] - yi)**2
        return Ey

    def position_constraint_y(y, *args):
        vertex = args[0]
        center = args[1]
        angle = args[2]
        J = args[3]
        yi = center[2] + y[0]
        for j in range(1, J + 1):
            yi += (y[2 * j - 1] * math.sin(j * angle) + y[2 * j] * math.cos(j * angle))
        yi -= vertex['position y']
        return yi

    def tangential_constraint_y(x, *args):
        vertex = args[0]
        angle = args[2]
        J = args[3]
        ti = 0
        for j in range(1, J+1):
            ti += (x[2 * j - 1] * j * math.cos(j * angle) - x[2 * j] * j * math.sin(j * angle))
        ti -= vertex['tangent y']
        return ti

    args1 = (previous, center, angles[0], J)
    args2 = (next, center, angles[-1], J)
    cons = ({'type': 'eq', 'fun': position_constraint_y, 'args': args1},
            {'type': 'eq', 'fun': position_constraint_y, 'args': args2},
            {'type': 'eq', 'fun': tangential_constraint_y, 'args': args1},
            {'type': 'eq', 'fun': tangential_constraint_y, 'args': args2}
            )
    y0 = np.zeros((2 * J + 1))  # initialize value
    res_y = minimize(fun_y, y0, args = (vertices, center, angles, J), method='SLSQP', constraints=cons)
    """
    print 'minimum value:'
    print res_y.fun
    print 'optimal resolution:'
    print res_y.x
    print 'is iteration sucessful:'
    print res_y.success
    print 'the reason for the termination of iteration:'
    print res_y.message
    """
    res_x = getCoefficients_for_end_composite_single(J, vertices, center, angles, previous, next, 0)
    return res_x.x, res_y.x


def getCoefficients_for_end_composite_single(J, Vertices, Center, Angles, previous, next, axis): # J is a scalar, not a vector
    I = len(Vertices)
    e = 1e-10  # a value that is very close to zero
    def fun_s(s, *args):  # function to be minimized
        Es = 0
        vertices = args[0]
        center = args[1]
        angles = args[2]
        JJ = args[3]
        for i in range(1, I - 1):
            si = center[axis] + s[0]
            for j in range(1, JJ + 1):
                si += (s[2 * j - 1] * math.cos(j * angles[i]) + s[2 * j] * math.sin(j * angles[i]))
            Es += (vertices[i][axis] - si)**2
        return Es

    def position_constraint_s(s, *args):
        vertex = args[0]
        center = args[1]
        angle = args[2]
        JJ = args[3]
        si = center[axis] + s[0]
        for j in range(1, JJ + 1):
            si += (s[2 * j - 1] * math.cos(j * angle) + s[2 * j] * math.sin(j * angle))
        if axis == 0:
            si -= vertex['position x']
        elif axis == 1:
            si -= vertex['position y']
        elif axis == 2:
            si -= vertex['position z']
        return si

    def tangential_constraint_s(s, *args):
        vertex = args[0]
        angle = args[2]
        JJ = args[3]
        ti = 0
        for j in range(1, JJ+1):
            ti += (s[2 * j - 1] * (-j) * math.sin(j * angle) + s[2 * j] * j * math.cos(j * angle))
        if axis == 0:
            ti -= vertex['tangent x']
        elif axis == 1:
            ti -= vertex['tangent y']
        elif axis == 2:
            ti -= vertex['tangent z']
        return ti

    args1 = (previous, Center, Angles[0], J)
    args2 = (next, Center, Angles[-1], J)
    cons = ({'type': 'eq', 'fun': position_constraint_s, 'args': args1},
            {'type': 'eq', 'fun': position_constraint_s, 'args': args2},
            {'type': 'eq', 'fun': tangential_constraint_s, 'args': args1},
            {'type': 'eq', 'fun': tangential_constraint_s, 'args': args2}
            )
    s0 = np.zeros((2 * J + 1))  # initialize value
    res_s = minimize(fun_s, s0, args=(Vertices, Center, Angles, J), method='SLSQP', constraints=cons)
    """
    print 'minimum value:'
    print res_s.fun
    print 'optimal resolution:'
    print res_s.x
    print 'is iteration successful:'
    print res_s.success
    print 'the reason for the termination of iteration:'
    print res_s.message
    """

    return res_s.x


def getCoefficients_for_end_composite_3D(J, vertices, center, angles, previous, next):
    coe = []
    for axis in [0, 1, 2]:
        res = getCoefficients_for_end_composite_single(J[axis], vertices, center, angles, previous, next, axis)
        coe.append(res)
    return coe


def getCoefficients_for_non_end_composite_single(J, vertices, center, angles, previous, axis): # J is a scalar, not a vector
    I = len(vertices)
    if len(center) != 3:
        raise ValueError("bad center, you input {}".format(center))
    if len(angles) != I:
        raise ValueError("the number of angles({}) does not equal to that of vertices{}".format(len(angles), I))
    if not axis in [0,1,2]:
        raise ValueError("axis should be one of 0,1,2")
    if isinstance(J, list):
        J = J[axis]
    elif isinstance(J, int):
        pass
    ConstArray = np.zeros(2 * J + 1)
    CoefficientMatrix = np.zeros(shape=(2 * J + 1, 2 * J + 1), dtype=float, order='C')  # row-major
    Pv0 = previous['position x']
    Tv0 = previous['tangent x']
    if axis == 1:
        Pv0 = previous['position y']
        Tv0 = previous['tangent y']
    elif axis == 2:
        Pv0 = previous['position z']
        Tv0 = previous['tangent z']
    v0 = angles[0]
    if v0 != 0.0 or v0 != math.pi:
        ConstArray[0] = Pv0 - center[axis] + Tv0 * math.cos(v0) / math.sin(v0)
        ConstArray[1] = -Tv0 / math.sin(v0)

        CoefficientMatrix[0, 0] = 1.
        CoefficientMatrix[0, 1] = 0.

        CoefficientMatrix[1, 0] = 0.
        CoefficientMatrix[1, 1] = 1.

        for j in range(2, J + 1):
            CoefficientMatrix[0, 2 * j - 1] = -j * math.sin(j * v0) * math.cos(v0) / math.sin(v0) + math.cos(j * v0)
            CoefficientMatrix[1, 2 * j - 1] = j * math.sin(j * v0) / math.sin(v0)

        for j in range(1, J + 1):
            CoefficientMatrix[0, 2 * j] = j * math.cos(j * v0) * math.cos(v0) / math.sin(v0) + math.sin(j * v0)
            CoefficientMatrix[1, 2 * j] = - j * math.cos(j * v0) / math.sin(v0)

        M = lambda jjj, v: jjj * math.sin(jjj * v0) * math.cos(v0) / math.sin(v0) - math.cos(jjj * v0) - jjj * math.cos(
            v) * math.sin(jjj * v0) / math.sin(v0) + math.cos(jjj * v)
        N = lambda jjj, v: -jjj * math.cos(jjj * v0) * math.cos(v0) / math.sin(v0) - math.sin(jjj * v0) + jjj * math.cos(
            jjj * v0) * math.cos(v) / math.sin(v0) + math.sin(jjj * v)

        for i in range(1, I):
            vi = angles[i]

            G = Pv0 + Tv0 / math.sin(v0) * (math.cos(v0) - math.cos(vi)) - vertices[i][axis]
            for k in range(2, 2 * J + 1):
                if k % 2 == 1:
                    m = M((k + 1) / 2, vi)
                    ConstArray[k] -= G * m
                    for jj in range(2, J + 1):
                        CoefficientMatrix[k, 2 * jj - 1] += M(jj, vi) * m
                    for jj in range(1, J + 1):
                        CoefficientMatrix[k, 2 * jj] += N(jj, vi) * m
                else:
                    n = N(k / 2, vi)
                    ConstArray[k] -= G * n
                    for jj in range(2, J + 1):
                        CoefficientMatrix[k, 2 * jj - 1] += M(jj, vi) * n
                    for jj in range(1, J + 1):
                        CoefficientMatrix[k, 2 * jj] += N(jj, vi) * n

    coe = np.linalg.solve(CoefficientMatrix, ConstArray)

    return coe


def getCoefficients_for_non_end_composite_3D(J, vertices, center, angles, previous):
    coe=[]
    for axis in range(3):
        coe.append(getCoefficients_for_non_end_composite_single(J[axis], vertices, center, angles, previous, axis))
    return coe


def getCoefficients_for_non_end_composite_2D(J, vertices, center, angles, previous):
    I = len(vertices)

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.zeros(shape=(2 * J + 1, 2 * J + 1), dtype=float, order='C')

    Pv0_y = previous['position y']
    Tv0_y = previous['tangent y']
    v0 = angles[0]

    if v0 != 0.0 or v0 != math.pi:
        bConstArray[0] = Pv0_y - center[2] - Tv0_y * math.sin(v0) / math.cos(v0)
        bConstArray[1] = Tv0_y / math.cos(v0)

        bCoefficientMatrix[0, 0] = 1.
        bCoefficientMatrix[0, 1] = 0.

        bCoefficientMatrix[1, 0] = 0.
        bCoefficientMatrix[1, 1] = 1.

        for j in range(2, J + 1):
            bCoefficientMatrix[0, 2 * j - 1] = - j * math.cos(j * v0) * math.sin(v0) / math.cos(v0) + math.sin(j * v0)
            bCoefficientMatrix[1, 2 * j - 1] = j * math.cos(j * v0) / math.cos(v0)

        for j in range(1, J + 1):
            bCoefficientMatrix[0, 2 * j] = j * math.sin(j * v0) * math.sin(v0) / math.cos(v0) + math.cos(j * v0)
            bCoefficientMatrix[1, 2 * j] = - j * math.sin(j * v0) / math.cos(v0)

        M_y = lambda j, v: j * math.cos(j * v0) * math.sin(v0) / math.cos(v0) - math.sin(j * v0) - j * math.cos(j * v0) * math.sin(v) / math.cos(v0) + math.sin(j * v)
        N_y = lambda j, v: -j * math.sin(j * v0) * math.sin(v0) / math.cos(v0) - math.cos(j * v0) + j * math.sin(j * v0) * math.sin(v) / math.cos(v0) + math.cos(j * v)

        for i in range(1, I):
            vi = angles[i]
            G_y = Pv0_y + Tv0_y / math.cos(v0) * (math.sin(vi) - math.sin(v0)) - vertices[i][2]
            for k in range(2, 2 * J + 1):
                if k % 2 == 1:
                    my = M_y((k + 1) / 2, vi)
                    bConstArray[k] -= G_y * my
                    for j in range(2, J + 1):
                        bCoefficientMatrix[k, 2 * j - 1] += M_y(j, vi) * my
                    for j in range(1, J + 1):
                        bCoefficientMatrix[k, 2 * j] += N_y(j, vi) * my
                else:
                    ny = N_y(k / 2, vi)
                    bConstArray[k] -= G_y * ny
                    for j in range(2, J + 1):
                        bCoefficientMatrix[k, 2 * j - 1] += M_y(j, vi) * ny
                    for j in range(1, J + 1):
                        bCoefficientMatrix[k, 2 * j] += N_y(j, vi) * ny

    b = np.linalg.solve(bCoefficientMatrix, bConstArray)

    a = getCoefficients_for_non_end_composite_single(J, vertices, center, angles, previous, axis=0)

    return a, b


def form_vertices_of_symmetry_ellipse_2D(vertices_first_half, center_first_half, angles_first_half, d_bar, a1, b1, a2, b2):
    I = len(vertices_first_half)
    numPt = len(d_bar)
    symmetry_ellipse_vertices = [[0 for i in range(3)] for j in range(numPt)]
    Ea = 0.0
    Em = 0.0
    d = [0.0] * numPt
    J1 = (len(a1) - 1) / 2
    J2 = (len(a2) - 1) / 2

    for i in range(I):
        new_first_half_vertex = [0.0] * 3
        new_second_half_vertex = [0.0] * 3

        new_first_half_vertex[0] = center_first_half[0] + a1[0]
        new_first_half_vertex[2] = center_first_half[2] + b1[0]
        v1 = angles_first_half[i]

        v2 = -angles_first_half[i]
        new_second_half_vertex[0] = -center_first_half[0] + a2[0]
        new_second_half_vertex[2] = center_first_half[2] + b2[0]

        for j in range(1, J1 + 1):
            new_first_half_vertex[0] += a1[2 * j - 1] * math.cos(j * v1) + a1[2 * j] * math.sin(j * v1)
            new_first_half_vertex[2] += b1[2 * j - 1] * math.sin(j * v1) + b1[2 * j] * math.cos(j * v1)
        for j in range(1, J2 + 1):
            new_second_half_vertex[0] += a2[2 * j - 1] * math.cos(j * v2) + a2[2 * j] * math.sin(j * v2)
            new_second_half_vertex[2] += b2[2 * j - 1] * math.sin(j * v2) + b2[2 * j] * math.cos(j * v2)

        symmetry_ellipse_vertices[i][0] = new_first_half_vertex[0]
        symmetry_ellipse_vertices[i][2] = new_first_half_vertex[2]
        di1 = math.sqrt((vertices_first_half[i][0] - new_first_half_vertex[0]) ** 2 + (vertices_first_half[i][2] - new_first_half_vertex[2]) ** 2)
        d[i] = di1
        Ea += (di1 / d_bar[i])
        if Em < di1 / d_bar[i]:
            Em = di1 / d_bar[i]

        if i != 0 and i != (I - 1):
            symmetry_ellipse_vertices[numPt - i][0] = new_second_half_vertex[0]
            symmetry_ellipse_vertices[numPt - i][2] = new_second_half_vertex[2]
            di2 = math.sqrt(pow(-vertices_first_half[i][0] - new_second_half_vertex[0], 2) + pow(vertices_first_half[i][2] - new_second_half_vertex[2], 2))
            d[numPt - i] = di2
            Ea += (di2 / d_bar[numPt - i])

            if Em < di2 / d_bar[numPt - i]:
                Em = di2 / d_bar[numPt - i]

    Ea = Ea / numPt

    return symmetry_ellipse_vertices, Ea, Em


def form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index):
    I = len(vertices)
    numPt = len(d_bar)
    fragment_vertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    J = (len(a) - 1) / 2
    tmp = form_vertices_of_fragment_single(a, center, angles, 0)

    for i in range(I):
        v_i = angles[i]
        y_i = center[2] + b[0]

        for j in range(1, J + 1):
            y_i += (b[2 * j - 1] * math.sin(j * v_i) + b[2 * j] * math.cos(j * v_i))

        x_i = tmp[i]
        fragment_vertices[i][0] = x_i
        fragment_vertices[i][2] = y_i

        d_i = math.sqrt((vertices[i][0] - x_i)**2 + (vertices[i][2] - y_i)**2)
        d.append(d_i)
        index = (i + start_index) % numPt
        Ea += (d_i / d_bar[index])
        if Em < d_i / d_bar[index]:
            Em = d_i / d_bar[index]

    Ea = Ea / I

    return fragment_vertices, Ea, Em


def form_vertices_of_fragment_3D(coe, vertices, center, angles, d_bar, start_index):
    I = len(vertices)
    numPt = len(d_bar)
    fragment_vertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    for axis in range(3):
        s = form_vertices_of_fragment_single(coe[axis], center, angles, axis)
        for i in range(I):
            fragment_vertices[i][axis] = s[i]

    for i in range(I):
        d_i = 0
        for axis in range(3):
            d_i = d_i + (vertices[i][axis] - fragment_vertices[i][axis])**2
        d_i = math.sqrt(d_i)
        d.append(d_i)
        index = (i + start_index) % numPt
        Ea += (d_i / d_bar[index])
        if Em < d_i / d_bar[index]:
            Em = d_i / d_bar[index]

    Ea = Ea / I

    return fragment_vertices, Ea, Em


# for both 2D and 3D
def form_vertices_of_fragment_single(coe, center, angles, axis): # coe is 1-dimensional
    I = len(angles)
    fragment_vertices = []

    if isinstance(coe[0], float):
        pass
    elif isinstance(coe[0], list):
        coe = coe[axis]
    J = (len(coe) - 1) / 2
    for i in range(I):
        v_i = angles[i]
        s_i = center[axis] + coe[0]

        for j in range(1, J + 1):
            s_i += (coe[2 * j - 1] * math.cos(j * v_i) + coe[2 * j] * math.sin(j * v_i))

        fragment_vertices.append(s_i)

    return fragment_vertices


# whole single curve, first segment, curve fragment
def findJ_2D(vertices,angles,d_bar,center,Ea_criteria,Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index):
    J = 0
    a3, b3 = func_getCoefficients(3, vertices, center, angles)
    v3, Ea3, Em3 = func_formGeneralizedEllipse(a3, b3, vertices, center, angles, d_bar,index)
    a10, b10 = func_getCoefficients(10, vertices, center, angles)
    v10, Ea10, Em10 = func_formGeneralizedEllipse(a10, b10, vertices, center, angles, d_bar,index)
    if Ea3 < Ea_criteria and Em3 < Em_criteria:
        J = find_smaller_J(vertices, angles, d_bar, center, 3, 10, Ea3, Ea10, Ea_criteria, Em_criteria, func_getCoefficients, func_formGeneralizedEllipse, index)
    elif Ea10 >= Ea_criteria or Em10 >= Em_criteria:
        J = find_bigger_J(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria, Em_criteria, func_getCoefficients, func_formGeneralizedEllipse, index)
    elif Ea3 >= Ea_criteria or Em3 >= Em_criteria:
        J = find_inbetween_J(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria,Em_criteria, func_getCoefficients, func_formGeneralizedEllipse, index)
    return J


def calculate_Ea_Em(before_vertices, after_vertices, d_bar, start_index, axis):
    Ea = 0
    Em = 0
    numPt = len(d_bar)
    for i in range(len(after_vertices)):
        d_i = 0
        N = 0
        if isinstance(after_vertices[0], float):
            d_i = abs(before_vertices[i][axis] - after_vertices[i])
            N = len(after_vertices)
        elif isinstance(after_vertices[0], om.MVector):
            d_i = abs(before_vertices[i][axis] - after_vertices[i][axis])
            N = len(after_vertices[0])
        index = (i + start_index) % numPt
        Ea += (d_i / d_bar[index])
        if Em < d_i / d_bar[index]:
            Em = d_i / d_bar[index]
    Ea = Ea / N
    return Ea, Em


# whole single curve, and first segment
def findJ_3D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, start_index):
    J = [3, 1, 3]
    coe = []

    for axis in [0, 1, 2]:# x,y,z axis
        Ea = 99999.9
        Em = 99999.9
        while Ea >= Ea_criteria/1.4 or Em >= Em_criteria:
            tmp = getCoefficients_single(J, vertices, center, angles, axis)
            fragment_vertices_i = form_vertices_of_fragment_single(tmp, center, angles, axis)
            Ea, Em = calculate_Ea_Em(vertices, fragment_vertices_i, d_bar, start_index, axis)
            J[axis] += 1

        coe.append(tmp)

    return J, coe


def find_inbetween_J(vertices,angles,d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria,Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index):
    # Linear interpolate to find J_small<J<J_big
    if J_small > J_big or J_small == J_big:
        raise ValueError('J_small({}) is no smaller than J_big({})'.format(J_small, J_big))

    if J_small == J_big - 1:
        return J_big

    J = 0
    if Ea_smallJ >= Ea_criteria:
        J = J_big - int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ))
    elif Em_smallJ >= Em_criteria:
        J = J_big - int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ))

    if J < J_small:
        raise ValueError('J({}) is smaller than J_small({})'.format(J, J_small))
    if J > J_big:
        raise ValueError('J({}) is bigger than J_big({})'.format(J, J_big))

    if J == J_small:
        J = J_small+1
    if J == J_big:
        J = J_big-1

    a, b = func_getCoefficients(J, vertices, center, angles)
    v, Ea, Em = func_formGeneralizedEllipse(a, b, vertices, center, angles, d_bar,index)
    if Ea < Ea_criteria and Em < Em_criteria: # J meets the criteria
        if J == J_small+1:
            return J
        else:
            return find_inbetween_J(vertices, angles, d_bar, center, J_small, J, Ea_smallJ, Em_smallJ, Ea, Em, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index)
    else:# J doesn't meet the criteria
        if J == J_big-1:
            return J_big
        else:
            return find_inbetween_J(vertices, angles, d_bar, center, J, J_big, Ea, Em, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index)


def find_bigger_J(vertices, angles, d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index):
    # Linear extrapolate to find bigger J>J_small
    if Ea_bigJ >= Ea_criteria:
        J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    elif Em_bigJ >= Em_criteria:
        J = int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ)) + J_small
    if J == J_big:
        J = J + 1
    if J < J_big:
        J = J_big + 1
    a, b = func_getCoefficients(J, vertices, center, angles)
    v, Ea, Em = func_formGeneralizedEllipse(a, b, vertices, center, angles, d_bar,index)
    if Ea >= Ea_criteria or Em >= Em_criteria:
        return find_bigger_J(vertices,angles,d_bar, center, J_big, J, Ea_bigJ, Em_bigJ, Ea, Em, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea < Ea_criteria and Em < Em_criteria and J > J_small:
            J -= 1
            a, b = func_getCoefficients(J, vertices, center, angles)
            v, Ea, Em = func_formGeneralizedEllipse(a, b, vertices, center, angles, d_bar,index)

        return J + 1


def find_smaller_J(vertices, angles, d_bar, center, J_small, J_big, Ea_smallJ, Ea_bigJ, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index):
    # Linear extrapolate to find smaller J<J_small<J_big,
    # The criteria is always Ea_criteria,
    # because both (Ea_smallJ,Em_smallJ) and (Ea_bigJ,Em_bigJ) meet criteria.
    # In this case, we will always use the average error Ea for the extrapolation since the average error is a global measurement.

    J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small

    if J<3:
        J = 3
        return J
    if J==J_small:
        return J_small
    a, b = func_getCoefficients(J, vertices, center, angles)
    v, Ea, Em = func_formGeneralizedEllipse(a, b, vertices, center, angles, d_bar,index)
    if Ea < Ea_criteria and Em < Em_criteria:
        return find_smaller_J(vertices, angles, d_bar, center, J, J_small, Ea, Ea_smallJ, Ea_criteria, Em_criteria,func_getCoefficients,func_formGeneralizedEllipse,index)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea >= Ea_criteria or Em >= Em_criteria and J < J_small:
            J += 1
            a, b = func_getCoefficients(J, vertices, center, angles,func_getCoefficients,func_formGeneralizedEllipse)

        return J


def findJ_for_non_end_composite_2D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, previous):
    J = 1
    start_index = previous['cut point index']
    a, b = getCoefficients_for_non_end_composite_2D(J, vertices, center, angles, previous)
    v, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index)

    while Ea >= Ea_criteria/1.5 or Em >= Em_criteria/1.5:
        J += 1
        a, b = getCoefficients_for_non_end_composite_2D(J, vertices, center, angles, previous)
        v, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index)

    return J


def findJ_for_non_end_composite_3D(vertices_one_segment, angles, d_bar, segment_center, Ea_criteria, Em_criteria, previous, func_getCoefficients_for_non_end_composite_single_3D):
    J = [1,1,1]
    coe = []
    start_index = previous['cut point index']
    for axis in [0, 1, 2]:  # x,y,z axis
        Ea = 99999.9
        Em = 0
        tmp = []
        for iterate in range(10):#while Ea >= (Ea_criteria/1.5) or Em >= (Em_criteria/1.5):
            tmp = func_getCoefficients_for_non_end_composite_single_3D(J, vertices_one_segment, segment_center, angles, previous, axis)
            fragment_vertices_i = form_vertices_of_fragment_single(tmp, segment_center, angles, axis)
            Ea, Em = calculate_Ea_Em(vertices_one_segment, fragment_vertices_i, d_bar, start_index, axis)
            if Ea <= (Ea_criteria/2.5) and Em <= (Em_criteria/2.5):
                break
            else:
                J[axis] += 1
        coe.append(tmp)

    return J, coe


def findJ_for_end_segment_2D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, previous, next):
    # for the end segment that links the first segment
    J = 3
    start_index = previous['cut point index']
    a, b = getCoefficients_for_end_composite_2D(3, vertices, center, angles, previous, next)
    v, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index)

    if Ea < Ea_criteria and Em < Em_criteria:
        while Ea < Ea_criteria and Em < Em_criteria and J > 1:
            J -= 1
            a, b = getCoefficients_for_end_composite_2D(J, vertices, center, angles, previous, next)
            v, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index)

        return J+1
    elif Ea >= Ea_criteria or Em >= Em_criteria:
        while Ea >= Ea_criteria or Em >= Em_criteria:
            J += 1
            a, b = getCoefficients_for_end_composite_2D(J, vertices, center, angles, previous, next)
            v, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices, center, angles, d_bar, start_index)

        return J


def findJ_for_end_segment_3D(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, previous, next):
    # for the end segment that links the first segment
    J = [1,1,1]
    start_index = previous['cut point index']
    for axis in [0,1,2]:
        Ea = 99999.9
        Em = 99999.9

        while Ea >= Ea_criteria/2.5 or Em >= Em_criteria/2.5:
            J[axis] += 1
            coe = getCoefficients_for_end_composite_single(J[axis], vertices, center, angles, previous, next, axis)
            v = form_vertices_of_fragment_single(coe, center, angles, axis)
            Ea, Em = calculate_Ea_Em(vertices, v, d_bar, start_index, axis)
    return J


def position_and_tangent_of_parametric_point_2D(a_adjacent, b_adjacent, angle):
    x_tan = 0
    y_tan = 0
    x = a_adjacent[0]
    y = b_adjacent[0]
    J = (len(a_adjacent)-1)/2
    for j in range(1, J+1):
        x_tan = x_tan + a_adjacent[2 * j] * j * math.cos(j * angle) - a_adjacent[2 * j - 1] * j * math.sin(j * angle)
        y_tan = y_tan + b_adjacent[2 * j - 1] * j * math.cos(j * angle) - b_adjacent[2 * j] * j * math.sin(j * angle)
        x += a_adjacent[2 * j - 1] * math.cos(j * angle) + a_adjacent[2 * j] * math.sin(j * angle)
        y += b_adjacent[2 * j - 1] * math.sin(j * angle) + b_adjacent[2 * j] * math.cos(j * angle)

    return x_tan, y_tan, x, y


def position_and_tangent_of_parametric_point_3D(coe_adjacent, angle):
    tan = [0, 0, 0]
    if len(coe_adjacent) != 3:
        raise ValueError("coe_adjacent should be a 3XN 2 dimentional list")
    if not isinstance(coe_adjacent[0][0], float):
        raise ValueError("coe_adjacent should be a 3XN 2 dimentional list of float")
    pos = [coe_adjacent[0][0], coe_adjacent[1][0], coe_adjacent[2][0]]
    J = [(len(coe_adjacent[0])-1)/2, (len(coe_adjacent[1])-1)/2, (len(coe_adjacent[2])-1)/2]
    for axis in range(3):
        for j in range(1, J[axis]+1):
            tan[axis] += coe_adjacent[axis][2 * j] * j * math.cos(j * angle) - coe_adjacent[axis][2 * j - 1] * j * math.sin(j * angle)
            pos[axis] += coe_adjacent[axis][2 * j - 1] * math.cos(j * angle) + coe_adjacent[axis][2 * j] * math.sin(j * angle)
    return tan, pos


def compisite_segment_with_one_end_shared_2D(index,vertices_matrix,angles_matrix,composite_a,composite_b,segment_center_list,
                                          d_bar,Ea_criteria, Em_criteria,cut_points,isClosed):
    x0_tan, y0_tan, x0, y0 = position_and_tangent_of_parametric_point_2D(composite_a[index - 1], composite_b[index - 1], angles_matrix[index - 1][-1])
    x0 += segment_center_list[index - 1][0]
    y0 += segment_center_list[index - 1][2]
    cut_pt_index = cut_points[index - 1]
    if isClosed == True:
        cut_pt_index = cut_points[index]
    previous = {'position x': x0, 'position y': y0, 'tangent x': x0_tan, 'tangent y': y0_tan,
                'cut point index': cut_pt_index}
    J = findJ_for_non_end_composite_2D(vertices_matrix[index], angles_matrix[index],d_bar, segment_center_list[index],
                                                  Ea_criteria, Em_criteria, previous)
    a, b = getCoefficients_for_non_end_composite_2D(J, vertices_matrix[index],segment_center_list[index],angles_matrix[index], previous)
    vertices, Ea, Em = form_vertices_of_fragment_2D(a, b, vertices_matrix[index],segment_center_list[index],angles_matrix[index],d_bar, cut_pt_index)

    return J,vertices,a,b,Ea,Em


def rebuild_curve(file_path, vertices, delete_points_list):
    # To rebuild the curve:
    # delete the unnecessary points
    # approximate the new list of points with trigonometric functions
    dir_path = cmds.workspace(fn=True) + '/data/'
    file_path = file_path.split('/')[-1]
    del vertices[delete_points_list.get(file_path)[2] + 1:delete_points_list.get(file_path)[3]]
    del vertices[delete_points_list.get(file_path)[0] + 1:delete_points_list.get(file_path)[1]]
    with open(dir_path + 'removed_' + file_path, "w") as f:
        for v in vertices:
            f.write('{} {} {}\n'.format(v[0], v[1], v[2]))
    # rebuild the curve for every cross-section
    center = getCenter_3D(vertices)
    angles = calculateAngle_3D(vertices, center)
    J = [20, 2, 20]
    coe = getCoefficients_3D(J, vertices, center, angles)
    # add points to replace the deleted points
    delta_angle = (angles[delete_points_list.get(file_path)[0] + 1] - angles[delete_points_list.get(file_path)[0]]) / (
                          delete_points_list.get(file_path)[1] - delete_points_list.get(file_path)[0])
    for index in range(delete_points_list.get(file_path)[0],
                       delete_points_list.get(file_path)[1] - 1):
        angles.insert(index + 1, delta_angle + angles[index])
    for index in range(delete_points_list.get(file_path)[2],
                       delete_points_list.get(file_path)[3] - 1):
        angles.insert(index + 1, delta_angle + angles[index])
    tmp_x = form_vertices_of_fragment_single(coe[0], center, angles, 0)
    tmp_y = form_vertices_of_fragment_single(coe[1], center, angles, 1)
    tmp_z = form_vertices_of_fragment_single(coe[2], center, angles, 2)
    new_vertices = []
    for i in range(len(tmp_x)):
        new_vertices.append(om.MVector(tmp_x[i], tmp_y[i], tmp_z[i]))
    vertices = new_vertices

    return vertices


if __name__ == "__main__":
    pass

