import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd
import quaternion
from mpl_toolkits.mplot3d import Axes3D

# 機体のスペックシートを読み込む
spec = pd.read_csv("spec.csv")
# 機体の推力シートを読む
thrust = np.loadtxt("G75.txt")
qua = quaternion.Quaternion()

T = 60.
N = 2000
dt = T / float(N)
t = np.empty(N)
t[0] = 0.
b = 1. / 6.
count = -1

# 機体の初期重量[kg]
m0 = float(spec['VALUE'][0])
# 機体の燃焼後重量[kg]
mb = float(spec['VALUE'][1])
# 先頭から見た初期重心の位置[m]
rg = float(spec['VALUE'][2])
# ボディチューブの長さ[m]
l0 = float(spec['VALUE'][3])
# エンジンの長さ[m]
le = float(spec['VALUE'][4])
# ノーズの長さ[m]
lc = float(spec['VALUE'][5])
m_cone = float(spec['VALUE'][6])
m_body = float(spec['VALUE'][7])

# 機体自体の初期重心[m](エンジンの部分を抜いた重心)
r = (m0 - mb) * (rg - le / 2.) / m0 + rg
# エンジン重心[m]
re = le / 2.

# ボディの内径[m]
ri = float(spec['VALUE'][8])
# ボディの外径[m]
r0 = float(spec['VALUE'][9])
# 空気抵抗係数
cd0 = float(spec['VALUE'][10])
# 圧力中心位置[m]
CP = float(spec['VALUE'][11])
# 重力加速度[m/s^2]
g = float(spec['VALUE'][12])
# パラシュート直径[m]
rp = float(spec['VALUE'][13])

# 燃焼時間
tt = thrust[:, 0]
# パラシュート直径
S = np.pi * r0 * r0 / 4.

fta = thrust[:, 1]
fta = np.r_[fta, 0.]
tn = np.r_[tt, T]
ft = interpolate.interp1d(tn, fta)

t1 = np.array([tt[0], tt[tt.size - 1], T])
m1 = np.array([m0, mb, mb])
m = interpolate.interp1d(t1, m1)

"""初期条件"""
q0 = np.array([np.cos(np.pi * 10. / 180.), np.sqrt(2) * np.sin(np.pi * 10. / 180.) / 2.,
               - np.sqrt(2) * np.sin(np.pi * 10. / 180.) / 2., 0.])
q = q0
p = np.empty([N, 3])
p[0] = np.array([0., 0., 0.])
v = np.empty([N, 3])
v[0] = np.array([0., 0., 0.])
a = np.empty([N, 3])
a[0] = np.array([0., 0., 0.])
vnorm = np.empty(N)
vnorm[0] = 0
anorm = np.empty(N)
anorm[0] = 0


def wind(h):
    wind_t = np.array([0., 0., 0.]) * pow((abs(h) / 10.), 1. / 7.)
    return wind_t


# パラシュート
<<<<<<< HEAD
Sp = np.pi * rp ** 2 / 4. - np.pi * rp ** 2 / 400
=======
Sp = np.pi * rp**2 / 4. - np.pi * rp**2 / 400
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
cp = 0.5

"""ωに関する設定"""
omega = np.empty([N, 3])
omega[0] = np.array([0., 0., 0.])
alpha = np.empty(N)
alpha[0] = 0
cd = np.empty(N)
cd[0] = cd0
kappa = np.empty(N)
kappa[0] = 1.293 * S * cd[0] / 2.

"""回転に関する関数"""


# 時間tでの機体の重心位置
def rcg(time):
    return (m(time) * r + (m(time) - mb) * (l0 + lc - re)) / (2. * m(time) - mb)


# ノーズの先端まわりの慣性モーメント（円錐と近似）
# http://ayapin.film.s.dendai.ac.jp/~matuda/Lecture/SRCS/phys11/inertiaofcone.pdfの平行軸なし
def cone():
    return (3. * m_cone * r0 ** 2 / 4.) / 20. + (3 * m_cone * lc ** 2) / 5.


# ボディの先端xy軸まわりの慣性モーメント(くり抜かれた円柱と近似)
# 重心x軸まわりの慣性モーメントhttps://www.orientalmotor.co.jp/tech/reference/sizing_motor03/
# そこから円錐の先頭中心に平行軸の定理で移動(最後まとめて移動させたいから)
def columnar():
    I_colum = m_body * ((r0 ** 2 + ri ** 2) / 4. + (lc ** 2) / 3.) / 4.
    dis = m_body * ((l0 / 2. + lc) ** 2 - (l0 / 2.) ** 2)
    return I_colum + dis


# エンジンの先端xy軸まわりの慣性モーメント（もう質点とみなしちゃう）[端からエンジン重心le/2.]
def inertia_e(time):
<<<<<<< HEAD
    return (m(time) - mb) * (lc + l0 - le / 2.) ** 2
=======
    return (m(time) - mb) * (lc+l0 - le/2.)**2
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf


# 最後平行軸の定理で全部ずらす(質量の合計はm(time))
def inertia_xy(time):
    return cone() + columnar() + inertia_e(time) + m(time) * rcg(time) * rcg(time)


# 円筒のz軸重心まわりの慣性モーメント
def inertia_z(time):
<<<<<<< HEAD
    return (r0 ** 2 + ri ** 2) * m(time) / 8.


# 慣性モーメントの時間微分を差分で
=======
    return (r0**2 + ri**2) * m(time) / 8.


>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
def inertia_xy_dot(time):
    return (inertia_xy(time + dt) - inertia_xy(time)) / dt


def inertia_z_dot(time):
    return (inertia_z(time + dt) - inertia_z(time)) / dt


<<<<<<< HEAD
# 圧力中心と重心の差(モーメント用)
=======
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
def length(time):
    return CP - rcg(time)


# 慣性モーメント設定
def I(time):
    MI = np.array([[inertia_xy(time), 0., 0.],
                   [0., inertia_xy(time), 0.],
                   [0., 0., inertia_z(time)]])
    return MI


<<<<<<< HEAD
# Iの逆行列の設定
=======
# 逆行列の設定
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
def Iinv(time):
    Mii = np.linalg.inv(I(time))
    return Mii


# 慣性モーメントの微分
def I_dot(time):
    Mid = np.array([[inertia_xy_dot(time), 0., 0.],
                    [0., inertia_xy_dot(time), 0.],
                    [0., 0., inertia_z_dot(time)]])
    return Mid


<<<<<<< HEAD
# 機体から見た抗力(ベクトルの向きは機体から見てるので-vb+wind)
=======
# 機体から見た抗力
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
def Fd(vb, a, h):
    drag = (- vb + wind(h)) * np.array([kappa[a] * np.linalg.norm(-vb + wind(h))])
    return drag


# 機体座標上の力のモーメント
def M(time, vb, qr, a, h):
    rv = np.array([0., 0., -length(time)])
    F = Fd(vb, a, h)
    Fl = qua.rotation(F, qr)
    Nmm = np.cross(rv, Fl)
    return Nmm


# 慣性モーメントと角速度の内積
def Iomega(time, a):
    return np.dot(I(time), omega[a])


# 慣性モーメント微分と角速度の内積
def Idomega(time, a):
    return np.dot(I_dot(time), omega[a])


# 角速度と「慣性モーメントと角速度の内積」の外積
<<<<<<< HEAD
def omegacross(time, a):
    return np.cross(omega[a], Iomega(time, a))


# θの二回微分=Frotの形にするために定義
def Frot(time, qr, vb, a, h):
    Frot_1 = np.dot(Iinv(time), M(time, vb, qr, a, h))
    Frot_2 = - np.dot(Iinv(time), Idomega(time, a))
    Frot_3 = - np.dot(Iinv(time), omegacross(time, a))
    return Frot_1 + Frot_2 + Frot_3


"""並進運動に関する関数"""


# 基軸座標系の推力を並進運動座標系に回した上で合力/mを考える
=======
def omegacross(a, time):
    return np.cross(omega[a], Iomega(time, a))


"""並進運動に関する関数"""


# 基軸座標系の推力を並進運動座標系に回す(抗力を回していない)
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
def Fs(time, qr, vb, a, h):
    W = np.array([0., 0., - m(time) * g])
    Ft = np.array([0., 0., ft(time)])
    Ftt = qua.rotation(Ft, qua.cquat(qr))
    Fa = Fd(vb, a, h)
<<<<<<< HEAD
=======
    # FD = qua.rotation(Fa, qua.cquat(qr))
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
    F = W + Ftt + Fa
    return F / m(time)


# 迎角
def angle(vb, qr):
    vr = (0., 0., 1.)
    vair = qua.rotation(vb, qr)
    theta = np.arccos(np.dot(vair, vr) / (np.linalg.norm(vair) * np.linalg.norm(vr)))
    return theta


# 計算
def runge_kutta(a, kg1, kg2, kg3, kg4):
    k1 = kg1 * np.array([b * dt])
    k2 = kg2 * np.array([2. * b * dt])
    k3 = kg3 * np.array([2. * b * dt])
    k4 = kg4 * np.array([b * dt])
    return a + k1 + k2 + k3 + k4


for i in range(N - 1):
    t[i + 1] = t[i] + dt
    """回転でωを求める"""
<<<<<<< HEAD
    ko1 = Frot(t[i], q, v[i], i, p[i, 2])
    kv1 = Fs(t[i], q, v[i], i, p[i, 2])
    kz1 = v[i]

    ko2 = Frot(t[i] + dt / 2., q, v[i] + kv1 * dt / 2., i, p[i, 2] + kz1[2] * dt / 2.)
    kv2 = Fs(t[i] + dt / 2., q, v[i] + kv1 * dt / 2., i, p[i, 2] + kz1[2] * dt / 2.)
    kz2 = v[i] + kv1 * dt / 2.

    ko3 = Frot(t[i] + dt / 2., q, v[i] + kv2 * dt / 2., i, p[i, 2] + kz2[2] * dt / 2.)
    kv3 = Fs(t[i] + dt / 2., q, v[i] + kv2 * dt / 2., i, p[i, 2] + kz2[2] * dt / 2.)
    kz3 = v[i] + kv2 * dt / 2.

    ko4 = Frot(t[i] + dt, q, v[i] + kv3 * dt, i, p[i, 2] + kz3[2] * dt / 2.)
=======
    ko1 = np.dot(Iinv(t[i]), M(t[i], v[i], q, i, p[i, 2])) - np.dot(Iinv(t[i]), Idomega(t[i], i)) - np.dot(Iinv(t[i]),
                                                                                                           omegacross(i,
                                                                                                                      t[
                                                                                                                          i]))
    kv1 = Fs(t[i], q, v[i], i, p[i, 2])
    kz1 = v[i]

    ko2 = np.dot(Iinv(t[i] + dt / 2.),
                 M(t[i] + dt / 2., v[i] + kv1 * dt / 2., q, i, p[i, 2] + kz1[2] * dt / 2.)) - np.dot(
        Iinv(t[i] + dt / 2.),Idomega(t[i] + dt / 2., i)) - np.dot(Iinv(t[i] + dt / 2.), omegacross(i, t[i] + dt / 2.))
    kv2 = Fs(t[i] + dt / 2., q, v[i] + kv1 * dt / 2., i, p[i, 2] + kz1[2] * dt / 2.)
    kz2 = v[i] + kv1 * dt / 2.

    ko3 = np.dot(Iinv(t[i] + dt / 2.),
                 M(t[i] + dt / 2., v[i] + kv2 * dt / 2., q, i, p[i, 2] + kz2[2] * dt / 2.)) - np.dot(
        Iinv(t[i] + dt / 2.),
        Idomega(t[i] + dt / 2.,
                i)) - np.dot(
        Iinv(t[i] + dt / 2.), omegacross(i, t[i] + dt / 2.))
    kv3 = Fs(t[i] + dt / 2., q, v[i] + kv2 * dt / 2., i, p[i, 2] + kz2[2] * dt / 2.)
    kz3 = v[i] + kv2 * dt / 2.

    ko4 = np.dot(Iinv(t[i] + dt), M(t[i] + dt, v[i] + kv3 * dt, q, i, p[i, 2] + kz3[2] * dt / 2.)) - np.dot(
        Iinv(t[i] + dt),
        Idomega(t[i] + dt, i)) - np.dot(Iinv(t[i] + dt), omegacross(i, t[i] + dt))
>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
    kv4 = Fs(t[i] + dt, q, v[i] + kv3 * dt, i, p[i, 2] + kz3[2] * dt / 2.)
    kz4 = v[i] + kv3 * dt

    omega[i + 1] = runge_kutta(omega[i], ko1, ko2, ko3, ko4)
    v[i + 1] = runge_kutta(v[i], kv1, kv2, kv3, kv4)
    p[i + 1] = runge_kutta(p[i], kz1, kz2, kz3, kz4)
    alpha[i + 1] = angle(v[i + 1], q)
    a[i] = kv1
    vnorm[i] = np.linalg.norm(v[i])
    anorm[i] = np.linalg.norm(a[i])
    print(p[i])

    if p[i + 1, 2] > p[i, 2]:
        cd[i + 1] = cd0 / abs(np.cos(alpha[i + 1]))

        if p[i + 1, 2] < 2.1:
            q = q0
            # ランチャー上を移動
            kappa[i + 1] = 1.293 * S * cd[i + 1] / 2.

        else:
            # クォータニオンを求める
            q += qua.qua_dot(omega[i], q) * dt
            kappa[i + 1] = 1.293 * S * cd[i + 1] / 2.

    else:
        # クォータニオンを求める
        q += qua.qua_dot(omega[i], q) * dt
        kappa[i + 1] = 1.293 * Sp * cp / 2.

    if p[i + 1, 2] < - 1.:
        count = i + 1
        break

<<<<<<< HEAD
=======

>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
plt.title("motion")
plt.xlabel("t[s]")
plt.ylim(0, max(p[0:count, 2] * 1.05))
plt.ylabel("z[m]")
plt.plot(t[0:count], p[0:count, 2])
plt.grid()
plt.show()

<<<<<<< HEAD
=======

>>>>>>> 750624a40987c73563fe4c26962a9afd5471ccbf
fig = plt.figure()
ax = Axes3D(fig)
ax.plot3D(p[0:count, 0], p[0:count, 1], p[0:count, 2])
plt.show()
