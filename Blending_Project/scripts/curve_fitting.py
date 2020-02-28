import math,sys, os
import operator
from PySide2 import QtCore, QtGui
import maya.api.OpenMaya as om

sys.path.append('/usr/lib64/python2.7/site-packages')
sys.path.append('./.local/lib/python2.7/site-packages')
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


def maya_useNewAPI():
    pass


class IllegalArgumentError(ValueError):
    pass


def getCenter(vertices):
    center = QtCore.QPointF(0.0, 0.0)
    for v in vertices:
        center += QtCore.QPointF(float(v[0]), float(v[2]))
    center = center / len(vertices)
    return center


def getCenter3(vertices):
    center = om.MVector(0.0, 0.0, 0.0)
    for v in vertices:
        center += v
    center = center / len(vertices)
    return center


def get_d_bar(vertices,center):
    d_bar = []
    for v in vertices:
        d_bar.append(math.sqrt((v[0] - center.x()) ** 2 + (v[2] - center.y()) ** 2))
    return d_bar


def get_d_bar3(vertices,center):
    d_bar = []
    for v in vertices:
        d_bar.append(math.sqrt((v[0] - center[0]) ** 2 + (v[1] - center[1]) ** 2 + (v[2] - center[2]) ** 2))
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

        if anglem > 2 * math.pi:
            anglem = anglem - 2 * math.pi

        angles.append(anglem)
    return angles


def calculateAngle3(vertices,center):
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
        dt = angles[(i + 1) % numPt] - angles[i - 1]
        if dt < 0:
            dt += math.pi*2
        tan_x = (vertices[(i + 1) % numPt][0] - vertices[i - 1][0]) / dt
        tan_y = (vertices[(i + 1) % numPt][2] - vertices[i - 1][2]) / dt
        tangents.append((tan_x, tan_y))
    return tangents


def calculateNormal(tangents):
    normals=[]
    for t in tangents:
        normal = QtGui.QVector2D(-t[1], t[0]).normalized()
        normals.append(normal)
    return normals

def calculateCurvature(tangents, angles):
    curvatures = []
    numPt = len(tangents)
    for i in range(0, numPt):
        dt = angles[(i + 1) % numPt] - angles[i - 1]
        if dt<0:
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
        dict[i]=abs(curvatures[i])

    sorted_d = sorted(dict.items(), key=operator.itemgetter(1))
    return sorted_d


def getCoefficients(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    if J<3:
        raise IllegalArgumentError('J must be bigger than 3, you input {0}'.format(J))
    I = len(vertices)
    aConstArray = np.zeros(2 * J + 1)
    aCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(2 * J + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

    for i in range(I):
        aCoefficientMatrix[0, i] = 1.
        bCoefficientMatrix[0, i] = 1.

    for i in range(I):  # for aCoefficientMatrix's column
        vi = angles[i]
        for j in range(1, J + 1):  # for aCoefficientMatrix's row
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


def getCoefficients3(J,vertices,center,angles):# abtain a[2j+1] and b[2j+1]
    I = len(vertices)
    aConstArray = np.zeros(2 * J + 1)
    aCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')  # row-major

    Jy = 2
    bConstArray = np.zeros(2 * Jy + 1)
    bCoefficientMatrix = np.ndarray(shape=(2 * Jy + 1, I), dtype=float, order='C')

    cConstArray = np.zeros(2 * J + 1)
    cCoefficientMatrix = np.ndarray(shape=(2 * J + 1, I), dtype=float, order='C')

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

            cCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            cCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            cConstArray[2 * j - 1] += (vertices[i][2] - center[2]) * math.cos(vi * j)
            cConstArray[2 * j] += (vertices[i][2] - center[2]) * math.sin(vi * j)
        for j in range(1, Jy + 1):
            bCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
            bCoefficientMatrix[2 * j, i] = math.sin(vi * j)

            bConstArray[2 * j - 1] += (vertices[i][1] - center[1]) * math.cos(vi * j)
            bConstArray[2 * j] += (vertices[i][1] - center[1]) * math.sin(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)
    C = np.dot(cCoefficientMatrix, cCoefficientMatrix.transpose())
    c = np.linalg.solve(C, cConstArray)

    return a, b, c


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


def formGeneralizedEllipse3(a, b, c, vertices, center, angles, d_bar):
    # CoefficientMatrix
    I = len(vertices)
    generalisedEllipseVertices = [[0 for i in range(3)] for j in range(I)]
    Ea = 0.0
    Em = 0.0
    d = []
    J = len(a) / 2
    Jy = len(b) / 2
    for i in range(I):
        generalisedEllipseVertices[i][0] = center[0] + a[0]
        generalisedEllipseVertices[i][1] = center[1] + b[0]
        generalisedEllipseVertices[i][2] = center[2] + c[0]
        v = angles[i]
        for j in range(1, J + 1):
            generalisedEllipseVertices[i][0] += a[2 * j - 1] * math.cos(j * v) + a[2 * j] * math.sin(j * v)
            generalisedEllipseVertices[i][2] += c[2 * j - 1] * math.cos(j * v) + c[2 * j] * math.sin(j * v)
        for j in range(1, Jy + 1):
            generalisedEllipseVertices[i][1] += b[2 * j - 1] * math.sin(j * v) + b[2 * j] * math.cos(j * v)

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

    center_fragment = QtCore.QPointF(0.0, 0.0)
    for v in vertices_fragment:
        center_fragment.setX(center_fragment.x() + v[0])
        center_fragment.setY(center_fragment.y() + v[2])
    center_fragment.setX(center_fragment.x() / len(vertices_fragment))
    center_fragment.setY(center_fragment.y() / len(vertices_fragment))

    return angles_fragment, vertices_fragment, center_fragment, start_index, end_index


def getCoefficients_for_first_generalized_elliptic_segment(vertices_first_segment, angles_first_segment, center_first_segment, J):
    I = len(vertices_first_segment)
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
                vi = angles_first_segment[i]
            except IndexError:
                raise
            else:
                aCoefficientMatrix[2 * j - 1, i] = math.cos(vi * j)
                aCoefficientMatrix[2 * j, i] = math.sin(vi * j)

                # aConstAtrray[0] and bConstAtrray[0] always equal to 0 by definition!
                aConstArray[2 * j - 1] += (vertices_first_segment[i][0] - center_first_segment.x()) * math.cos(
                    vi * j)
                aConstArray[2 * j] += (vertices_first_segment[i][0] - center_first_segment.x()) * math.sin(
                    vi * j)

                bCoefficientMatrix[2 * j - 1, i] = math.sin(vi * j)
                bCoefficientMatrix[2 * j, i] = math.cos(vi * j)

                bConstArray[2 * j - 1] += (vertices_first_segment[i][2] - center_first_segment.y()) * math.sin(vi * j)
                bConstArray[2 * j] += (vertices_first_segment[i][2] - center_first_segment.y()) * math.cos(vi * j)

    A = np.dot(aCoefficientMatrix, aCoefficientMatrix.transpose())
    a = np.linalg.solve(A, aConstArray)
    B = np.dot(bCoefficientMatrix, bCoefficientMatrix.transpose())
    b = np.linalg.solve(B, bConstArray)

    return a, b


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


def getCoefficients_for_second_half_of_composite_segment(J, vertices_second_half, angles_second_half, center_first_half, center_second_half, a_first_half, b_first_half):
    I = len(vertices_second_half)
    J1 = (len(a_first_half)-1) / 2
    aConstArray = np.zeros(I)
    aCoefficientMatrix = np.ndarray(shape=(2 * J - 3, I), dtype=float, order='C')  # row-major

    bConstArray = np.zeros(I)
    bCoefficientMatrix = np.ndarray(shape=(2 * J - 3, I), dtype=float, order='C')
    for i in range(len(vertices_second_half)):
        vi = angles_second_half[i]
        Cai = vertices_second_half[i][0] - (center_first_half.x() + a_first_half[0]) * math.cos(2 * vi) - center_second_half.x() * (1 - math.cos(2 * vi))
        Cbi = vertices_second_half[i][2] - (center_first_half.y() + b_first_half[0]) * math.cos(2 * vi) - center_second_half.y() * (1 - math.cos(2 * vi))
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
    a3 = center_first_half.x() - center_second_half.x() + a_first_half[0] - a[0]
    b4 = center_first_half.y() - center_second_half.y() + b_first_half[0] - b[0]
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


def form_vertices_of_symmetry_ellipse(vertices_first_half, center_first_half, angles_first_half, d_bar, a1, b1, a2, b2):
    I = len(vertices_first_half)
    numPt = len(d_bar)
    symmetry_ellipse_vertices = [[0 for i in range(2)] for j in range(numPt)]
    Ea = 0.0
    Em = 0.0
    d = [0.0] * numPt
    J1 = (len(a1) - 1) / 2
    J2 = (len(a2) - 1) / 2

    for i in range(I):
        new_first_half_vertex = [0.0] * 2
        new_second_half_vertex = [0.0] * 2

        new_first_half_vertex[0] = center_first_half.x() + a1[0]
        new_first_half_vertex[1] = center_first_half.y() + b1[0]
        v1 = angles_first_half[i]

        v2 = -angles_first_half[i]
        new_second_half_vertex[0] = -center_first_half.x() + a2[0]
        new_second_half_vertex[1] = center_first_half.y() + b2[0]

        for j in range(1, J1 + 1):
            new_first_half_vertex[0] += a1[2 * j - 1] * math.cos(j * v1) + a1[2 * j] * math.sin(j * v1)
            new_first_half_vertex[1] += b1[2 * j - 1] * math.sin(j * v1) + b1[2 * j] * math.cos(j * v1)
        for j in range(1, J2 + 1):
            new_second_half_vertex[0] += a2[2 * j - 1] * math.cos(j * v2) + a2[2 * j] * math.sin(j * v2)
            new_second_half_vertex[1] += b2[2 * j - 1] * math.sin(j * v2) + b2[2 * j] * math.cos(j * v2)

        symmetry_ellipse_vertices[i][0] = new_first_half_vertex[0]
        symmetry_ellipse_vertices[i][1] = new_first_half_vertex[1]
        di1 = math.sqrt((vertices_first_half[i][0] - new_first_half_vertex[0]) ** 2 + (vertices_first_half[i][2] - new_first_half_vertex[1]) ** 2)
        d[i] = di1
        Ea += (di1 / d_bar[i])
        if Em < di1 / d_bar[i]:
            Em = di1 / d_bar[i]

        if i != 0 and i != (I - 1):
            symmetry_ellipse_vertices[numPt - i][0] = new_second_half_vertex[0]
            symmetry_ellipse_vertices[numPt - i][1] = new_second_half_vertex[1]
            di2 = math.sqrt(pow(-vertices_first_half[i][0] - new_second_half_vertex[0], 2) + pow(vertices_first_half[i][2] - new_second_half_vertex[1], 2))
            d[numPt - i] = di2
            Ea += (di2 / d_bar[numPt - i])

            if Em < di2 / d_bar[numPt - i]:
                Em = di2 / d_bar[numPt - i]

    Ea = Ea / numPt

    return symmetry_ellipse_vertices, Ea, Em


def form_vertices_of_composite_ellipse(vertices_first_half, vertices_second_half, center_first_half, center_second_half, angles_first_half, angles_second_half, d_bar, a1, b1, a2, b2):
    I = len(vertices_first_half)
    numPt = len(d_bar)
    composite_ellipse_vertices = [[0 for i in range(2)] for j in range(numPt)]
    Ea = 0.0
    Em = 0.0
    d = [0.0] * numPt
    J1 = (len(a1) - 1) / 2
    J2 = (len(a2) - 1) / 2

    for i in range(I):
        new_first_half_vertex = [0.0] * 2
        new_second_half_vertex = [0.0] * 2

        new_first_half_vertex[0] = center_first_half.x() + a1[0]
        new_first_half_vertex[1] = center_first_half.y() + b1[0]
        v1 = angles_first_half[i]

        v2 = angles_second_half[i]
        new_second_half_vertex[0] = center_second_half.x() + a2[0]
        new_second_half_vertex[1] = center_second_half.y() + b2[0]

        for j in range(1, J1 + 1):
            new_first_half_vertex[0] += a1[2 * j - 1] * math.cos(j * v1) + a1[2 * j] * math.sin(j * v1)
            new_first_half_vertex[1] += b1[2 * j - 1] * math.sin(j * v1) + b1[2 * j] * math.cos(j * v1)
        for j in range(1, J2 + 1):
            new_second_half_vertex[0] += a2[2 * j - 1] * math.cos(j * v2) + a2[2 * j] * math.sin(j * v2)
            new_second_half_vertex[1] += b2[2 * j - 1] * math.sin(j * v2) + b2[2 * j] * math.cos(j * v2)

        composite_ellipse_vertices[i][0] = new_first_half_vertex[0]
        composite_ellipse_vertices[i][1] = new_first_half_vertex[1]
        di1 = math.sqrt((vertices_first_half[i][0] - new_first_half_vertex[0]) ** 2 + (vertices_first_half[i][2] - new_first_half_vertex[1]) ** 2)
        d[i] = di1
        Ea += (di1 / d_bar[i])
        if Em < di1 / d_bar[i]:
            Em = di1 / d_bar[i]

        if i != 0 and i != (I - 1):
            composite_ellipse_vertices[I -1 + i][0] = new_second_half_vertex[0]
            composite_ellipse_vertices[I -1 + i][1] = new_second_half_vertex[1]
            di2 = math.sqrt(pow(vertices_second_half[i][0] - new_second_half_vertex[0], 2) + pow(vertices_second_half[i][2] - new_second_half_vertex[1], 2))
            d[numPt - i] = di2
            Ea += (di2 / d_bar[I - 1 + i])

            if Em < di2 / d_bar[I - 1 + i]:
                Em = di2 / d_bar[I - 1 + i]

    Ea = Ea / numPt

    return composite_ellipse_vertices, Ea, Em


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


def form_vertices_of_fragment(a, b, vertices, center, angles, d_bar, start_index):
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

    Ea = Ea / I

    return fragment_vertices, Ea, Em


def findJ(vertices,angles,d_bar,center,Ea_criteria,Em_criteria):
    J=0
    a3,b3=getCoefficients(3, vertices, center, angles)
    v3,Ea3,Em3=formGeneralizedEllipse(a3, b3, vertices, center, angles, d_bar)
    a10,b10=getCoefficients(10, vertices, center, angles)
    v10,Ea10,Em10=formGeneralizedEllipse(a10, b10, vertices, center, angles, d_bar)
    if Ea3<Ea_criteria and Em3<Em_criteria:
        J=find_smaller_J(vertices, angles, d_bar, center, 3, 10, Ea3, Ea10, Ea_criteria, Em_criteria)
    elif Ea10>=Ea_criteria or Em10>=Em_criteria:
        J=find_bigger_J(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria, Em_criteria)
    elif Ea3>=Ea_criteria or Em3>=Em_criteria:
        J=find_inbetween_J(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria,Em_criteria)
    return J


def find_inbetween_J(vertices,angles,d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria,Em_criteria):
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

    a, b = getCoefficients(J, vertices, center, angles)
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea < Ea_criteria and Em < Em_criteria: # J meets the criteria
        if J == J_small+1:
            return J
        else:
            return find_inbetween_J(vertices, angles, d_bar, center, J_small, J, Ea_smallJ, Em_smallJ, Ea, Em, Ea_criteria, Em_criteria)
    else:# J doesn't meet the criteria
        if J == J_big-1:
            return J_big
        else:
            return find_inbetween_J(vertices, angles, d_bar, center, J, J_big, Ea, Em, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria)


def find_bigger_J(vertices, angles, d_bar, center,J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria):
    # Linear extrapolate to find bigger J>J_small
    if Ea_bigJ >= Ea_criteria:
        J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    elif Em_bigJ >= Em_criteria:
        J = int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ)) + J_small
    if J == J_big:
        J = J + 1
    if J < J_big:
        J = J_big + 1
    a, b = getCoefficients(J, vertices, center, angles)
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea >= Ea_criteria or Em >= Em_criteria:
        return find_bigger_J(vertices,angles,d_bar, center, J_big, J, Ea_bigJ, Em_bigJ, Ea, Em, Ea_criteria, Em_criteria)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea < Ea_criteria and Em < Em_criteria and J > J_small:
            J -= 1
            a, b = getCoefficients(J, vertices, center, angles)
            v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles,d_bar)

        return J + 1


def find_smaller_J(vertices, angles, d_bar, center, J_small, J_big, Ea_smallJ, Ea_bigJ, Ea_criteria, Em_criteria):
    # Linear extrapolate to find smaller J<J_small<J_big,
    # The criteria is always Ea_criteria,
    # because both (Ea_smallJ,Em_smallJ) and (Ea_bigJ,Em_bigJ) meet criteria.
    # In this case, we will always use the average error Ea for the extrapolation since the average error is a global measurement.

    J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    if J<3:
        J = 3
        return J
    a, b = getCoefficients(J, vertices, center, angles)
    v, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)
    if Ea < Ea_criteria and Em < Em_criteria:
        return fingSmallerJ(vertices, angles, d_bar, center, J, J_small, Ea, Ea_smallJ, Ea_criteria, Em_criteria)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea >= Ea_criteria or Em >= Em_criteria and J < J_small:
            J += 1
            a, b = getCoefficients(J, vertices, center, angles)

        return J


def findJ_for_end_segment(vertices, angles, d_bar, center, Ea_criteria, Em_criteria, previous, next):
    J = 5
    start_index = previous['cut point index']
    a3, b3 = getCoefficients_for_end_composite(3, vertices, center, angles, previous, next)
    v3, Ea3, Em3 = form_vertices_of_fragment(a3, b3, vertices, center, angles, d_bar, start_index)
    a10, b10 = getCoefficients_for_end_composite(10, vertices, center, angles, previous, next)
    v10, Ea10, Em10 = form_vertices_of_fragment(a10, b10, vertices, center, angles, d_bar, start_index)
    if Ea3 < Ea_criteria and Em3 < Em_criteria:
        J = find_smaller_J_for_end_composite(vertices, angles, d_bar, center, 3, 10, Ea3, Ea10, Ea_criteria, Em_criteria, previous, next)
    elif Ea10 >= Ea_criteria or Em10 >= Em_criteria:
        J = find_bigger_J_for_end_composite(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria, Em_criteria, previous, next)
    elif Ea3 >= Ea_criteria or Em3 >= Em_criteria:
        J = find_inbetween_J_for_end_composite(vertices, angles, d_bar, center, 3, 10, Ea3, Em3, Ea10, Em10, Ea_criteria, Em_criteria, previous, next)
    return J


def find_inbetween_J_for_end_composite(vertices, angles, d_bar, center, J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ,
                                                   Ea_criteria, Em_criteria, previous, next):
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
        J = J_small + 1
    if J == J_big:
        J = J_big - 1

    start_index = previous['cut point index']
    a, b = getCoefficients_for_end_composite(J, vertices, center, angles, previous, next)
    v, Ea, Em = form_vertices_of_fragment(a, b, vertices, center, angles, d_bar, start_index)
    if Ea < Ea_criteria and Em < Em_criteria:  # J meets the criteria
        if J == J_small + 1:
            return J
        else:
            return find_inbetween_J_for_end_composite(vertices, angles, d_bar, center, J_small, J, Ea_smallJ, Em_smallJ, Ea, Em,
                                    Ea_criteria, Em_criteria, previous, next)
    else:  # J doesn't meet the criteria
        if J == J_big - 1:
            return J_big
        else:
            return find_inbetween_J_for_end_composite(vertices, angles, d_bar, center, J, J_big, Ea, Em, Ea_bigJ, Em_bigJ, Ea_criteria,
                                    Em_criteria, previous, next)


def find_bigger_J_for_end_composite(vertices, angles, d_bar, center, J_small, J_big, Ea_smallJ, Em_smallJ, Ea_bigJ, Em_bigJ, Ea_criteria, Em_criteria, previous, next):
    # Linear extrapolate to find bigger J>J_small
    print 'find bigger J. J_small is {} and J_big is {}'.format(J_small, J_big)
    start_index = previous['cut point index']
    J = J_big
    if Ea_bigJ >= Ea_criteria:
        J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    elif Em_bigJ >= Em_criteria:
        J = int((Em_criteria - Em_smallJ) * (J_big - J_small) / (Em_bigJ - Em_smallJ)) + J_small
    if J <= J_big:
        print 'weird, J({}) is smaller than J_big({})'.format(J, J_big)
        J = J_big + 1
        return J

    a, b = getCoefficients_for_end_composite(J, vertices, center, angles, previous, next)
    v, Ea, Em = form_vertices_of_fragment(a, b, vertices, center, angles, d_bar, start_index)

    if Ea >= Ea_criteria or Em >= Em_criteria:
        find_bigger_J_for_end_composite(vertices, angles, d_bar, center, J_big, J, Ea_bigJ, Em_bigJ, Ea, Em, Ea_criteria, Em_criteria, previous, next)

    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea < Ea_criteria and Em < Em_criteria and J > J_small:
            J -= 1
            a, b = getCoefficients_for_end_composite(J, vertices, center, angles, previous, next)
            v, Ea, Em = form_vertices_of_fragment(a, b, vertices, center, angles, d_bar, start_index)

        return J + 1


def find_smaller_J_for_end_composite(vertices, angles, d_bar, center, J_small, J_big, Ea_smallJ, Ea_bigJ, Ea_criteria, Em_criteria, previous, next):
    # Linear extrapolate to find smaller J<J_small<J_big,
    # The criteria is always Ea_criteria,
    # because both (Ea_smallJ,Em_smallJ) and (Ea_bigJ,Em_bigJ) meet criteria.
    # In this case, we will always use the average error Ea for the extrapolation since the average error is a global measurement.

    J = int((Ea_criteria - Ea_smallJ) * (J_big - J_small) / (Ea_bigJ - Ea_smallJ)) + J_small
    if J < 3:
        J = 3
        return J

    start_index = previous['cut point index']
    a, b = getCoefficients_for_end_composite(J, vertices, center, angles, previous, next)
    v, Ea, Em = form_vertices_of_fragment(a, b, vertices, center, angles, d_bar, start_index)
    if Ea < Ea_criteria and Em < Em_criteria:
        return fing_smallerJ_for_end_composite(vertices, angles, d_bar, center, J, J_small, Ea, Ea_smallJ, Ea_criteria, Em_criteria, previous, next)
    else:
        # we are close to the solution, hence, a while function will suffice
        while Ea >= Ea_criteria or Em >= Em_criteria and J < J_small:
            J += 1
            a, b = getCoefficients_for_end_composite(J, vertices, center, angles, previous, next)

        return J
    return J


if __name__ == "__main__":
    """    
    file_paths = [
        'Source_Chest_cross_section_u_at_0_percentage_worldspace.dat',
        'Source_Chest_cross_section_u_at_10_percentage_worldspace.dat',
        'Source_Chest_cross_section_u_at_22_percentage_worldspace.dat'
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
        'Source_Head_cross_section_u_at_0_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_5_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_10_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_20_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_30_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_40_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_50_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_60_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_70_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_80_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_90_percentage_worldspace.dat',
        'Source_Head_cross_section_u_at_100_percentage_worldspace.dat'
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
        'Source_LeftArm_cross_section_u_at_35_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_40_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_50_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_55_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_60_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_70_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_80_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_90_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_95_percentage_worldspace.dat',
        'Source_LeftArm_cross_section_u_at_100_percentage_worldspace.dat'
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
    """
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

        file_name = file_path.replace('Source_', '')
        directory = dirPath + file_name.split('_cross_section_')[0]

        try:
            os.stat(directory)
        except:
            os.mkdir(directory)

        save_file_path = directory + '/' + file_name

        J = 6

        if 'worldspace' not in file_path:
            center = getCenter(vertices)
            d_bar = get_d_bar(vertices, center)
            angles = calculateAngle(vertices, center)
            a, b = getCoefficients(J, vertices, center, angles)
            generalisedEllipseVertices, Ea, Em = formGeneralizedEllipse(a, b, vertices, center, angles, d_bar)

            with open(save_file_path, "w+") as f:
                f.write('range:')
                f.write('0-360 \n')
                f.write('center: {} {}\n'.format(center.x(), center.y()))
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
            center = getCenter3(vertices)
            d_bar = get_d_bar3(vertices, center)
            angles = calculateAngle3(vertices, center)
            a, b, c = getCoefficients3(J, vertices, center, angles)
            generalisedEllipseVertices, Ea, Em = formGeneralizedEllipse3(a, b, c, vertices, center, angles, d_bar)

            with open(save_file_path, "w+") as f:
                f.write('range:')
                f.write('0-360 \n')
                f.write('center: {} {} {}\n'.format(center[0], center[1], center[2]))
                f.write('a: ')
                for i in a:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('b: ')
                for i in b:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('c: ')
                for i in c:
                    f.write(str(i) + ' ')
                f.write('\n')
                f.write('angles: ')
                for i in angles:
                   f.write(str(i)+' ')
                f.write('\n')
        """
        image_name = file_name.split('.dat')[0] + '_generalised_ellipse.png'
        image_dir = cmds.workspace(fn=True)+'/images/'+file_name.split('_cross_section_')[0]+ '/'
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
        save_image_path = image_dir + image_name
        
        file=QtCore.QFile(file_path)
        try:
            file.open(QtCore.QIODevice.WriteOnly):
            pixmap=QtGui.QPixmap((500,400))
            self.canvas.render(pixmap)
            pixmap.save(file,"PNG")
            file.close()
        """        


