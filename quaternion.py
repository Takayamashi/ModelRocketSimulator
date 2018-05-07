import numpy as np


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

    """クォータニオンの時間微分"""

    def qua_dot(self, omega, q):
        omega = self.quaternion(omega)
        return - 0.5 * self.qcross(omega, q)
