import numpy as np
from copy import deepcopy
from numpy.linalg import norm, det
from numpy import sin, cos, dot

class IK_Math:
    @staticmethod
    def Vmat(ang):
        ang = (3.1415926536 - ang)
        return np.array([[cos(ang), -sin(ang), 0, 0], [sin(ang), cos(ang), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    @staticmethod
    def Tmat(ang):
        ang = (-ang)
        return np.array([[1, 0, 0, 0], [0, cos(ang), -sin(ang), 0], [0, sin(ang), cos(ang), 0], [0, 0, 0, 1]])

    @staticmethod
    def Dmat(dis):
        return np.array([[1, 0, 0, dis], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    @staticmethod
    def xyz_from_frame(mat):
        return np.array([mat[0][3], mat[1][3], mat[2][3]])

    @staticmethod
    def dir_from_frame(mat):
        return np.array([mat[0][0], mat[1][0], mat[2][0]])

    @staticmethod
    def gettorsion(points):
        # 4 points are given
        fr1_side = points[0] - points[1]
        fr1_mid = points[2] - points[1]
        fr2_mid = -fr1_mid
        fr2_side = points[3] - points[2]

        fr1_side -= dot(fr1_side, fr1_mid) * fr1_mid / dot(fr1_mid, fr1_mid)
        fr2_side -= dot(fr2_side, fr2_mid) * fr2_mid / dot(fr2_mid, fr2_mid)

        for i in [fr1_side,fr2_side]:
            i /= norm(i)
        ang = np.arccos(np.dot(fr1_side, fr2_side))
        if det(np.array([fr1_side, fr1_mid, fr2_side])) < 0:
            ang = -ang
        return ang

    @staticmethod
    def gettorsion_vec(a, b, c):
        # 3 vectors are given
        fr1_side = a
        fr1_mid = b
        fr2_mid = -b
        fr2_side = c

        fr1_side -= dot(fr1_side, fr1_mid) * fr1_mid / dot(fr1_mid, fr1_mid)
        fr2_side -= dot(fr2_side, fr2_mid) * fr2_mid / dot(fr2_mid, fr2_mid)

        for i in [fr1_side, fr2_side]:
            i /= norm(i)
        ang = np.arccos(np.dot(fr1_side, fr2_side))
        if det(np.array([fr1_side, fr1_mid, fr2_side])) < 0:
            ang = -ang
        return ang

    @staticmethod
    def getvalangle(points):
        prevvec = points[0]
        curvec = points[1]
        nextvec = points[2]
        dir1 = prevvec - curvec
        dir2 = nextvec - curvec
        dir1 /= norm(dir1)
        dir2 /= norm(dir2)
        angle = np.arccos(np.dot(dir1, dir2))
        return angle

    @staticmethod
    def gs_rand(x, y):
        x /= norm(x)
        y -= np.dot(x, y) * x
        y /= norm(y)
        z = np.array([1.0, 1.0, 1.0])
        z -= np.dot(x, z) * x
        z -= np.dot(y, z) * y
        z /= norm(z)
        if det(np.array([x, y, z])) < 0:
            z = -z
        return z
