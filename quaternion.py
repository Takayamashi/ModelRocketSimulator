import numpy as np
"""
クォータニオンをベクトルq = [q0, q1, q2, q3]で定義する．
共役クォータニオンq* = [q0, -q1, -q2, -q3]
クォータニオン積を行列表記すると，
qp = [[q0, -q1, -q2, -q3], [p0,
      [q1,  q0, -q3,  q2],  p1,
      [q2,  q3,  q0, -q1],  p2,
      [q3, -q2,  q1,  q0]]  p3]
ある座標系のベクトルr_0をq_01によって座標系を回転させると，
r_1 = q_01* r_1 q_01
"""


class Quaternion:
    """クォータニオンに関する関数の定義"""

    # 三次元ベクトルをクォータニオンに変換
    def quaternion(self, u):
        q = np.r_[0., u]
        return q

    """共役クォータニオン"""

    def cquat(self, q):
        qt = np.array([q[0], -q[1], -q[2], -q[3]])
        return qt

    "クォータニオン積"

    def qcross(self, q1, q2):
        q1 = np.array([[q1[0], -q1[1], -q1[2], -q1[3]],
                       [q1[1], q1[0], -q1[3], q1[2]],
                       [q1[2], q1[3], q1[0], -q1[1]],
                       [q1[3], -q1[2], q1[1], q1[0]]])
        return np.dot(q1, q2)

    # ベクトルrをクォータニオンqで回す
    def rotation(self, r, q):
        r = self.quaternion(r)
        coq = self.cquat(q)
        ra = self.qcross(coq, self.qcross(r, q))
        return np.array([ra[1], ra[2], ra[3]])

    # ベクトルrをクォータニオンq*で回す
    def rotation_w_r(self, r, q):
        r = self.quaternion(r)
        coq = self.cquat(q)
        ra = self.qcross(q, self.qcross(r, coq))
        return np.array([ra[1], ra[2], ra[3]])


    """クォータニオンの時間微分"""

    def qua_dot(self, omega, q):
        omega = self.quaternion(omega)
        return - 0.5 * self.qcross(omega, q)
