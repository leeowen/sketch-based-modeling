import math,sys, os
from PySide2 import QtCore

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy as np


def maya_useNewAPI():
    pass


def getCenter(vertices):
    center = QtCore.QPointF(0.0, 0.0)
    for v in vertices:
        center += QtCore.QPointF(float(v[0]), float(v[2]))
    center = center / len(vertices)
    return center


def get_d_bar(vertices,center):
    d_bar = []
    for v in vertices:
        d_bar.append(math.sqrt((v[0] - center.x()) ** 2 + (v[2] - center.y()) ** 2))
    return d_bar


def calculateAngle(vertices,center):
    angles = []
    for v in vertices:
        anglem = (v[2] - center.y()) / math.sqrt(
            (v[2] - center.y()) ** 2 + (v[0] - center.x()) ** 2)
        anglem = math.acos(anglem)

        if v[0] > center.x() and v[2] > center.y():
            pass
        elif v[0] > center.x() and v[2] < center.y():
            pass
        elif v[0] < center.x() and v[2] > center.y():
            anglem = 2 * math.pi - anglem
        elif v[0] < center.x() and v[2] < center.y():
            anglem = 2 * math.pi - anglem

        angles.append(anglem)
    return angles


def calculateArcLength(vertices):
    totalArcLength = 0.0
    arcLength = []
    numPt = len(vertices)
    for i in range(0, numPt):
        arcL = math.sqrt(pow(vertices[(i + 1) % numPt][0] - vertices[i][0], 2) + pow(
            vertices[(i + 1) % numPt][2] - vertices[i][2], 2))
        arcLength.append(arcL)
        totalArcLength += arcLength
    return totalArcLength, arcLength


def calculateTangent(vertices, angles):
    tangents = []
    numPt = len(vertices)
    for i in range(0, numPt):
        # primes refer to the derivatives dx/dt, dy/dt with respect to the angle t
        tan_x = (vertices[(i + 1) % numPt][0] - vertices[i - 1][0]) / (
                    angles[(i + 1) % numPt] - angles[i - 1])
        tan_y = (vertices[(i + 1) % numPt][2] - vertices[i - 1][2]) / (
                    angles[(i + 1) % numPt] - angles[i - 1])
        tangents.append((tan_x, tan_y))
    return tangents


def calculateCurvature(tangents, angles):
    curvatures = []
    numPt = len(tangents)
    for i in range(0, numPt):
        d2xdt2 = (tangents[(i + 1) % numPt][0] - tangents[i][0]) / (
                    angles[(i + 1) % numPt] - angles[i - 1])
        d2ydt2 = (tangenst[(i + 1) % numPt][1] - tangents[i][1]) / (
                    angles[(i + 1) % numPt] - angles[i - 1])
        k = (tangents[i][0] * d2ydt2 - tangents[i][1] * d2xdt2) / pow(
            pow(tangents[i][0], 2) + pow(tangents[i][1], 2), 3.0 / 2.0)
        curvatures.append(k)
    return curvatures


def getCoefficients(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    I = len(vertices)
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
            aConstArray[2 * j - 1] += (vertices[i][0] - center.x()) * math.cos(vi * j)
            aConstArray[2 * j] += (vertices[i][0] - center.x()) * math.sin(vi * j)

            bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
            bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][2] - center.y()) * math.sin(vi * j)
            bConstArray[2 * j] += (vertices[i][2] - center.y()) * math.cos(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    return a, b


def formGeneralizedEllipse(a, b, vertices, center, angles, d_bar):
    # CoefficientMatrix
    I = len(vertices)
    generalisedEllipseVertices = [[0 for i in range(2)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    J = len(a) / 2
    for i in range(I):
        generalisedEllipseVertices[i][0] = center.x() + a[0]
        generalisedEllipseVertices[i][1] = center.y() + b[0]
        v = angles[i]
        for j in range(1, J + 1):
            generalisedEllipseVertices[i][0] += a[2 * j - 1] * math.cos(j * v) + a[2 * j] * math.sin(j * v)
            generalisedEllipseVertices[i][1] += b[2 * j - 1] * math.sin(j * v) + b[2 * j] * math.cos(j * v)

        di = math.sqrt((vertices[i][0] - generalisedEllipseVertices[i][0]) ** 2 + (
                    vertices[i][2] - generalisedEllipseVertices[i][1]) ** 2)
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

    center_fragment = QtCore.QPointF(0.0, 0.0)
    for v in vertices_fragment:
        center_fragment.setX(center_fragment.x() + v[0])
        center_fragment.setY(center_fragment.y() + v[2])
    center_fragment.setX(center_fragment.x() / len(vertices_fragment))
    center_fragment.setY(center_fragment.y() / len(vertices_fragment))

    return angles_fragment, vertices_fragment, center_fragment, start_index, end_index


def coefficients_solver_for_first_half_of_segmented_ellipse(vertices_first_half, angles_first_half, center_first_half, J):
    I = len(vertices_first_half)
    aConstArray = np.zeros(2 * J + 1)
    aCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

    for i in range(I):
        aCoefficientMatrix[0, i] = 1.
        bCoefficientMatrix[0, i] = 1.

    for i in range(I):  # for aCoefficientMatrix's column
        for j in range(1, J + 1):  # for aCoefficientMatrix's row
            try:
                vi = angles_first_half[i]
            except IndexError:
                raise
            else:
                aCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
                aCoefficientMatrix[2 * j, i] = math.sin(vi * j)

                # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                aConstArray[2 * j - 1] += (vertices_first_half[i][0] - center_first_half.x()) * math.cos(
                    vi * j)
                aConstArray[2 * j] += (vertices_first_half[i][0] - center_first_half.x()) * math.sin(
                    vi * j)

                bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
                bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

                bConstArray[2 * j - 1] += (vertices_first_half[i][2] - center_first_half.y()) * math.sin(vi * j)
                bConstArray[2 * j] += (vertices_first_half[i][2] - center_first_half.y()) * math.cos(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    return a, b


def coefficients_solver_for_second_half_of_segmented_ellipse(vertices_second_half, angles_second_half, center_first_half, center_second_half, J, a_first_half, b_first_half):
    I = len(vertices_second_half)

    Ca = np.zeros((I, 1), order='C')  # order='C' means row-major
    Cb = np.zeros((I, 1), order='C')

    Ma = np.zeros((2 * J - 3, I), order='C')
    Mb = np.zeros((2 * J - 3, I), order='C')

    for i in range(I):
        v_i = angles_second_half[i]
        tmp = 1 - math.cos(2 * v_i)
        Ma[0][i] = tmp
        Mb[0][i] = tmp
        x_i = vertices_second_half[i][0]
        y_i = vertices_second_half[i][2]
        Ca[i][0] = x_i - (center_first_half.x() + a_first_half[0]) * math.cos(
            2 * v_i) - center_second_half.x() * (1 - math.cos(2 * v_i))
        Cb[i][0] = y_i - (center_first_half.y() + b_first_half[0]) * math.cos(
            2 * v_i) - center_second_half.y() * (1 - math.cos(2 * v_i))
        for j in range(1, J + 1):
            D1 = (1 - pow(-1, j)) * math.cos(v_i) + (1 + pow(-1, j)) * math.cos(2 * v_i)
            D2 = (1 - pow(-1, j)) * math.sin(v_i) + 1 / 2.0 * (1 + pow(-1, j)) * math.sin(2 * v_i)
            Ca[i][0] -= 1 / 2.0 * (D1 * a_first_half[2 * j - 1] + j * D2 * a_first_half[2 * j])
            Cb[i][0] -= 1 / 2.0 * (j * D2 * b_first_half[2 * j - 1] + D1 * b_first_half[2 * j])
        for j in range(3, J + 1):
            D1 = (1 - pow(-1, j)) * math.cos(v_i) + (1 + pow(-1, j)) * math.cos(2 * v_i)
            tmp = math.cos(j * v_i) - 1 / 2.0 * D1
            Ma[j * 2 - 5][i] = tmp
            Mb[j * 2 - 4][i] = tmp
            D2 = (1 - pow(-1, j)) * math.sin(v_i) + 1 / 2.0 * (1 + pow(-1, j)) * math.sin(2 * v_i)
            tmp = math.sin(j * v_i) - j / 2.0 * D2
            Ma[j * 2 - 4][i] = tmp
            Mb[j * 2 - 5][i] = tmp

    A = Ma.dot(Ma.transpose())
    aConstArray = Ma.dot(Ca)
    a_tmp = np.linalg.solve(A, aConstArray)
    B = np.dot(Mb, Mb.transpose())
    bConstArray = np.dot(Mb, Cb)
    b_tmp = np.linalg.solve(B, bConstArray)

    # get a1,a2,a3,a4 and b1,b2,b3,b4
    a1 = 0
    a2 = 0
    a3 = center_first_half.x() - center_second_half.x() + a_first_half[0] - a_tmp[0]
    a4 = 0
    b1 = 0
    b2 = 0
    b3 = 0
    b4 = center_first_half.y() - center_second_half.y() + b_first_half[0] - b_tmp[0]

    for j in range(1, 3):
        a1 += 0.5 * (1 - pow(-1, j)) * a_first_half[2 * j - 1]
        a2 += 0.5 * (1 - pow(-1, j)) * j * a_first_half[2 * j]
        a3 += 0.5 * (1 + pow(-1, j)) * a_first_half[2 * j - 1]
        a4 += 0.25 * (1 + pow(-1, j)) * j * a_first_half[2 * j - 1]
        b1 += 0.5 * (1 - pow(-1, j)) * b_first_half[2 * j - 1]
        b2 += 0.5 * (1 - pow(-1, j)) * j * b_first_half[2 * j]
        b3 += 0.25 * (1 + pow(-1, j)) * j * b_first_half[2 * j - 1]
        b4 += 0.5 * (1 + pow(-1, j)) * b_first_half[2 * j - 1]
    for j in range(3, J + 1):
        a1 += 0.5 * (1 - pow(-1, j)) * (a_first_half[2 * j - 1] - a_tmp[2 * j - 5])
        a2 += 0.5 * (1 - pow(-1, j)) * j * (a_first_half[2 * j] - a_tmp[2 * j - 4])
        a3 += 0.5 * (1 + pow(-1, j)) * (a_first_half[2 * j - 1] - a_tmp[2 * j - 5])
        a4 += 0.25 * (1 + pow(-1, j)) * j * (a_first_half[2 * j] - a_tmp[2 * j - 4])
        b1 += 0.5 * (1 - pow(-1, j)) * j * (b_first_half[2 * j - 1] - b_tmp[2 * j - 5])
        b2 += 0.5 * (1 - pow(-1, j)) * (b_first_half[2 * j] - b_tmp[2 * j - 4])
        b3 += 0.25 * (1 + pow(-1, j)) * j * (b_first_half[2 * j - 1] - a_tmp[2 * j - 5])
        b4 += 0.5 * (1 + pow(-1, j)) * (b_first_half[2 * j] - b_tmp[2 * j - 4])

    a = np.zeros(2 * J + 1)
    b = np.zeros(2 * J + 1)
    a[0] = a_tmp[0]
    a[1] = a1
    a[2] = a2
    a[3] = a3
    a[4] = a4
    b[0] = b_tmp[0]
    b[1] = b1
    b[2] = b2
    b[3] = b3
    b[4] = b4
    for j in range(3, J + 1):
        a[2 * j - 1] = a_tmp[2 * j - 5]
        a[2 * j] = a_tmp[2 * j - 4]
        b[2 * j - 1] = b_tmp[2 * j - 5]
        b[2 * j] = b_tmp[2 * j - 4]

    return a, b


def form_vertices_of_segmented_ellipse(vertices_first_half, vertices_second_half, center_first_half,center_second_half,
                                       angles_first_half, angles_second_half, d_bar, a1, b1, a2, b2):
    I = len(vertices_first_half)
    numPt = len(d_bar)
    segmented_ellipse_vertices = [[0 for i in range(2)] for j in range(numPt)]
    Ea = 0.0
    Em = 0.0
    d = [0.0] * numPt
    J = (len(a1) - 1) / 2

    for i in range(I):
        new_first_half_vertex = [0.0] * 2
        new_second_half_vertex = [0.0] * 2

        new_first_half_vertex[0] = center_first_half.x() + a1[0]
        new_first_half_vertex[1] = center_first_half.y() + b1[0]
        v1 = angles_first_half[i]

        v2 = angles_second_half[i]
        new_second_half_vertex[0] = (center_first_half.x() + a1[0]) * math.cos(2 * v2) + (center_second_half.x() + a2[0]) * (1 - math.cos(2 * v2))
        new_second_half_vertex[1] = (center_first_half.y() + b1[0]) * math.cos(2 * v2) + (center_second_half.y() + b2[0]) * (1 - math.cos(2 * v2))

        for j in range(1, J + 1):
            new_first_half_vertex[0] += a1[2 * j - 1] * math.cos(j * v1) + a1[2 * j] * math.sin(j * v1)
            new_first_half_vertex[1] += b1[2 * j - 1] * math.sin(j * v1) + b1[2 * j] * math.cos(j * v1)
            D1 = (1 - pow(-1, j)) * math.cos(v2) + (1 + pow(-1, j)) * math.cos(2 * v2)
            D2 = (1 - pow(-1, j)) * math.sin(v2) + 0.5 * (1 + pow(-1, j)) * math.sin(2 * v2)
            new_second_half_vertex[0] += 0.5 * (D1 * a1[2 * j - 1] + j * D2 * a1[2 * j])
            new_second_half_vertex[1] += 0.5 * (j * D2 * b1[2 * j - 1] + D1 * b1[2 * j])
        for j in range(3, J + 1):
            D1 = (1 - pow(-1, j)) * math.cos(v2) + (1 + pow(-1, j)) * math.cos(2 * v2)
            D2 = (1 - pow(-1, j)) * math.sin(v2) + 0.5 * (1 + pow(-1, j)) * math.sin(2 * v2)
            new_second_half_vertex[0] += (math.cos(j * v2) - 0.5 * D1) * a2[2 * j - 1] + (
                        math.sin(j * v2) - j / 2.0 * D2) * a2[2 * j]
            new_second_half_vertex[1] += (math.sin(j * v2) - j / 2.0 * D2) * b2[2 * j - 1] + (
                        math.cos(j * v2) - D1 / 2.0) * b2[2 * j]

        segmented_ellipse_vertices[i][0] = new_first_half_vertex[0]
        segmented_ellipse_vertices[i][1] = new_first_half_vertex[1]
        di1 = math.sqrt((vertices_first_half[i][0] - new_first_half_vertex[0]) ** 2 + (vertices_first_half[i][2] - new_first_half_vertex[1]) ** 2)
        d[i] = di1
        Ea += (di1 / d_bar[i])
        if Em < di1 / d_bar[i]:
            Em = di1 / d_bar[i]

        if i != 0 and i != (I - 1):
            segmented_ellipse_vertices[i + I - 1][0] = new_second_half_vertex[0]
            segmented_ellipse_vertices[i + I - 1][1] = new_second_half_vertex[1]
            di2 = math.sqrt(pow(vertices_second_half[i][0] - new_second_half_vertex[0], 2) + pow(vertices_second_half[i][2] - new_second_half_vertex[1], 2))
            d[i + I - 1] = di2
            Ea += (di2 / d_bar[i + I - 1])

            if Em < di2 / d_bar[i + I - 1]:
                Em = di2 / d_bar[i + I - 1]

    Ea = Ea / numPt

    return segmented_ellipse_vertices, Ea, Em


def coefficients_solver_for_fragmented_ellipse( J, angles, vertices, center):
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
            aConstArray[2 * j - 1] += (vertices[i][0] - center.x()) * math.cos(vi * j)
            aConstArray[2 * j] += (vertices[i][0] - center.x()) * math.sin(vi * j)

            bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
            bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][2] - center.y()) * math.sin(vi * j)
            bConstArray[2 * j] += (vertices[i][2] - center.y()) * math.cos(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    return a, b


def form_vertices_of_fragment(a, b, vertices, angles, center, d_bar, start_index):
    I = len(vertices)
    numPt = len(d_bar)
    fragment_vertices = [[0 for i in range(2)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    J = (len(a) - 1) / 2

    for i in range(I):
        v_i = angles[i]
        x_i = center.x() + a[0]
        y_i = center.y() + b[0]

        for j in range(1, J + 1):
            x_i += a[2 * j - 1] * math.cos(j * v_i) + a[2 * j] * math.sin(j * v_i)
            y_i += b[2 * j - 1] * math.sin(j * v_i) + b[2 * j] * math.cos(j * v_i)

        fragment_vertices[i][0] = x_i
        fragment_vertices[i][1] = y_i
        d_i = math.sqrt(pow(vertices[i][0] - x_i, 2) + pow(vertices[i][2] - y_i, 2))
        d.append(d_i)
        index = (i + start_index) % numPt
        Ea += (d_i / d_bar[index])
        if Em < d_i / d_bar[index]:
            Em = d_i / d_bar[index]

    Ea = Ea / numPt

    return fragment_vertices, Ea, Em


def findJ(vertices,angles,d_bar,center,Ea_criteria,Em_criteria):
    J=0
    a3,b3=getCoefficients(3,vertices,center,angles)
    v3,Ea3,Em3=formGeneralizedEllipse(a3, b3, vertices, center, angles, d_bar)
    a10,b10=getCoefficients(10,vertices,center,angles)
    v10,Ea10,Em10=formGeneralizedEllipse(a10, b10, vertices, center, angles, d_bar)
    if Ea3<Ea_criteria and Em3<Em_criteria:
        J=find_smaller_J(vertices,angles,d_bar,center,3,10, Ea3, Ea10,Ea_criteria,Em_criteria)
    elif Ea10>=Ea_criteria or Em10>=Em_criteria:
        J=find_bigger_J(vertices,angles,d_bar, center,3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria, Em_criteria)
    elif Ea3>=Ea_criteria or Em3>=Em_criteria:
        J=find_inbetween_J(vertices,angles,d_bar, center,3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria,Em_criteria)
    return J


def find_inbetween_J(vertices,angles,d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ,Ea_bigJ, Em_bigJ, Ea_criteria,Em_criteria):
    # Linear interpolate to find J_small<J<J_big
    if J_small > J_big:
        raise ValueError('J_small({}) is bigger than J_big({})'.format(J_small, J_big))
    """
    if Ea_smallJ < Ea_criteria and Ea_bigJ >= Ea_criteria:
        tmp = Ea_smallJ
        Ea_smallJ = Ea_bigJ
        Ea_bigJ = tmp

    if Em_smallJ < Em_criteria and Em_bigJ >= Em_criteria:
        tmp = Em_smallJ
        Em_smallJ = Em_bigJ
        Em_bigJ = tmp
    """

    if J_small == J_big - 1:
        return J_big

    J = 0
    if Ea_smallJ >= Ea_criteria:
        J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    elif Em_smallJ >= Em_criteria:
        J = int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ)) + J_small

    if J == J_small:
        J += 1

    a, b = getCoefficients(J, vertices, center, angles)
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea < Ea_criteria and Em < Em_criteria:
        if J == J_small + 1:
            return J
        else:
            return find_inbetween_J(vertices,angles,d_bar, center, J_small, J, Ea_smallJ, Ea, Em_smallJ, Em, Ea_criteria, Em_criteria)
    elif Ea >= Ea_criteria or Em >= Em_criteria:
        if J == J_big - 1:
            return J_big
        else:
            return find_inbetween_J(vertices,angles,d_bar, center,J, J_big, Ea, Ea_bigJ, Em, Em_bigJ, Ea_criteria, Em_criteria)


def find_bigger_J(vertices,angles,d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria):
    # Linear extrapolate to find bigger J>J_small
    if Ea_bigJ >= Ea_criteria:
        J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    elif Em_bigJ >= Em_criteria:
        J = int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ)) + J_small
    if J > J_big:
        a, b = getCoefficients(J, vertices, center, angles)
    else:
        return J_big
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea >= Ea_criteria or Em >= Em_criteria:
        return find_bigger_J(vertices,angles,d_bar, center,J, J_big, Ea, Em, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea < Ea_criteria and Em < Em_criteria and J > J_small:
            J -= 1
            a, b = getCoefficients(J, vertices, center, angles)
            v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles,d_bar)

        return J + 1


def find_smaller_J(vertices,angles,d_bar,center,J_small, J_big, Ea_smallJ, Ea_bigJ,Ea_criteria,Em_criteria):
    # Linear extrapolate to find smaller J<J_small<J_big,
    # The criteria is always Ea_criteria,
    # because both (Ea_smallJ,Em_smallJ) and (Ea_bigJ,Em_bigJ) meet criteria.
    # In this case, we will always use the average error Ea for the extrapolation since the average error is a global measurement.

    J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    a, b = getCoefficients(J, vertices, center, angles)
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea < Ea_criteria and Em < Em_criteria:
        return fingSmallerJ(J, J_small, Ea, Ea_smallJ)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea >= Ea_criteria or Em >= Em_criteria and J < J_small:
            J += 1
            a, b = getCoefficients(J, vertices, center, angles)

        return J

