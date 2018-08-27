from pyexcel_xls import get_data
import json
from math import pi, sqrt, atan, degrees, cos, sin, atan2
from data import Data
import airfoils as af
import time
import copy
import threading


m0 = 0
R = 0.0
r0 = 0.0
b = 0.0
RPM = 0.0
wR = 0.0
Sigma = 0.0
z = 0
N = 30
Fi7 = 8
delk = 0
delv = 0
rho0 = 1.225
khi = 0.96
Iv = 1.1
Cvt = 1.02
r = []
_r_ = []
dr = []
dFi = []
dPsi = 1


def main():
    global twist
    global threads
    global start_time
    global blade_mass
    global deltaT
    threads = 3 #multiprocessing.cpu_count()
    pnts = N/threads
    data = get_data("test.xls")
    hdata = json.dumps(data['HeliData'])
    twist = json.dumps(data['Twist'])
    blade_mass = json.dumps(data['Mass'])
    blade_calc(blade_mass)
    get_helidata(hdata)
    r.append(r0)
    _r_.append(r0 / R)
    dr.append(0)
    i = 1
    af.read()

    while i < N:
        r.append(r[i - 1] + (R - r0) / (N - 1))
        _r_.append(r[i]/R)
        dr.append(r[i] - r[i - 1])
        i = i + 1
    twst = get_twist(twist)
    twists(twst)
    deltaT = (dPsi/360) * 2 * pi / (wR / R)
    #print(deltaT)
    start_time = time.time()
    dt = init(Fi7)
    print("---  Init time %s seconds ---" % (time.time() - start_time))
    Q = []
    psi = dPsi
    start_time = time.time()
    Q.append(copy.deepcopy(dt))
    for i in range(int(360 / (dPsi))):
        calc(Fi7, Q[i], psi)
        Q.append(copy.deepcopy(Q[i]))
        psi = psi + dPsi

    print("---  Calc time %s seconds ---" % (time.time() - start_time))

def thrd(Q, s_psi, n):
    psi = dPsi + s_psi
    for i in range(int(120/dPsi)):
        print("z=%d\t" % (n), end=' ')
        calc(Fi7, Q, psi)
        psi = psi + dPsi


def get_helidata(hdata):
    js = json.loads(hdata)
    rt = dict()
    global r0
    global m0
    global R
    global b
    global wR
    global K
    global z
    global Sigma
    global lhh
    global dPsiy
    for q in js:
        rt.update({q[0] : q[1]})
    m0 = rt['m0']
    R = rt['R']
    r0 = rt['r0']
    b = rt['b']
    z = rt['z']
    K = rt['k']
    lhh = rt['lhh']
    RPM = rt['RPM']
    dPsiy = rt['dPsi']
    wR = 202  # RPM*pi*R/30  #После отладки исправить
    Sigma = b*z/(pi*R)


def get_twist(twist):
    js = json.loads(twist)
    rt = dict()
    for q in js[1:]:
        rt.update({q[1] : q[2]})
    return rt


def aprox(Arg, X, Y):
    i = 0
    j = 0
    while i < len(X):
        if X[j] <= Arg:
            j = j + 1
        i = i + 1
    x = Arg
    x1 = X[j - 1]
    x2 = X[j]
    y1 = Y[j - 1]
    y2 = Y[j]
    return (x - x1) * (y2 - y1) / (x2 - x1) + y1


def calc(Fi7, dt, psi = 0, delv = 0, delk = 0, z = 1):
    start_time1 = time.time()
    dt.Vyt = dt.Vy / wR
    dt.Fi7 = Fi7
    rho = rho0 * (20000 - dt.H) / (20000 + dt.H)
    dpsi = psi - dt.psi
    i = 0
    o_beta = dt.beta
    o_dbeta = dt.dbeta
    o_d2beta = dt.d2beta
    dbeta = o_dbeta+o_d2beta*dpsi
    dt.dbeta = dbeta
    beta = o_beta + dbeta*dpsi
    dt.beta = beta
    dt.psi = psi
    dt.T = 0
    dt.Mki = 0
    dt.Mkp = 0
    dt.ST = 0
    dt.V_y = dt.Vy/wR
    dST = []
    while i < N:
        dt.G0[i] = (dt.G0[i]+dt.G1[i])/2
        dt.u1[i] = dt.G0[i]/_r_[i]
        dt.v1[i] = -dt.V_y/2 + sqrt((dt.V_y/2)**2+abs(dt.G0[i]))*dt.G0[i]/abs(dt.G0[i])
        dt.V1[i] = dt.v1[i] + dt.V_y + (r[i] - lhh)*dt.dbeta/R
        dt.U1[i] = _r_[i] - dt.u1[i]
        dt.W1[i] = sqrt(dt.V1[i]**2+dt.U1[i]**2)
        dt.f[i] = atan(dt.V1[i]/dt.U1[i])
        dt.fd[i] = degrees(dt.f[i])
        dt.Fi[i] = dt.Fi7+dFi[i]-K*degrees(dt.beta)
        dt.Alpha[i] = dt.Fi[i] - dt.fd[i]
        dt.M[i] = dt.W1[i]*wR/340
        dt.Cy[i] = af.findCy(dt.Alpha[i], dt.M[i])
        dt.Cx[i] = af.findCx(dt.Alpha[i], dt.M[i])
        dt.dy[i] = 0.5*dt.Cy[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dx[i] = 0.5*dt.Cx[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dt[i] = (dt.dy[i]*cos(dt.f[i])-dt.dx[i]*sin(dt.f[i]))*cos(dt.beta)
        dt.dq[i] = dt.dx[i]*cos(dt.f[i])+dt.dy[i]*sin(dt.f[i])
        dt.dmki[i] = dt.dy[i]*sin(dt.f[i])*r[i]*Iv/khi
        dt.dmkp[i] = dt.dx[i]*cos(dt.f[i])*r[i]*Cvt
        dt.G1[i] = Sigma*dt.Cy[i]*dt.W1[i]/8
        dST.append(dt.dt[i] * (r[i] - lhh))
        i = i + 1

    i = 1

    while i < N:
        dt.dY[i] = (dt.dy[i]+dt.dy[i-1])*dr[i-1]/2
        dt.ST = dt.ST + (dST[i] + dST[i - 1]) * dr[i - 1] / 2
        dt.dX[i] = (dt.dx[i]+dt.dx[i-1])*dr[i-1]/2
        dt.dT[i] = (dt.dt[i]+dt.dt[i-1])*dr[i-1]/2
        dt.dQ[i] = (dt.dq[i]+dt.dq[i-1])*dr[i-1]/2
        dt.dMki[i] = (dt.dmki[i]+dt.dmki[i-1])*dr[i-1]/2
        dt.dMkp[i] = (dt.dmkp[i] + dt.dmkp[i - 1]) * dr[i - 1] / 2
        dt.T = dt.T + dt.dT[i]
        dt.Mki = dt.Mki + dt.dMki[i]
        dt.Mkp = dt.Mkp + dt.dMkp[i]
        i = i + 1
    i = 0
    A1 = (Ihh*cos(dt.beta)-lhh*Shh)*(wR/R)**2*sin(dt.beta)
    dt.d2beta = (dt.ST - A1-9.81*Shh)/(Ihh*(wR/R)**2)
    #print(dt.d2beta)

    dt.T = dt.T  * khi
    dt.Mki = dt.Mki
    dt.Mkp = dt.Mkp
    dt.Mk = dt.Mki + dt.Mkp
    print("psi = %f\tbeta = %f\tT = %f\tMk = %f\t%f\t time = %f" % (psi, degrees(dt.beta) ,dt.T, dt.Mk, dt.d2beta, time.time() - start_time))
    return dt


def init(Fi7):
    dt = Data(N, Fi7)
    dt.Vyt = dt.Vy / wR
    rho = rho0*(20000-dt.H)/(20000+dt.H)
    i = 0
    dST = []
    dt.V_y = dt.Vy/wR
    while i < N:
        dt.u1[i] = dt.G0[i]/_r_[i]
        dt.v1[i] = -dt.V_y/2 + sqrt((dt.V_y/2)**2+abs(dt.G0[i]))*dt.G0[i]/abs(dt.G0[i])
        dt.V1[i] = dt.v1[i] + dt.V_y
        dt.U1[i] = _r_[i] - dt.u1[i]
        dt.W1[i] = sqrt(dt.V1[i]**2+dt.U1[i]**2)
        dt.f[i] = atan(dt.V1[i]/dt.U1[i])
        dt.fd[i] = degrees(dt.f[i])
        dt.Fi[i] = dt.Fi7+dFi[i]-K*dt.beta
        dt.Alpha[i] = dt.Fi[i] - dt.fd[i]
        dt.M[i] = dt.W1[i]*wR/340
        dt.Cy[i] = af.findCy(dt.Alpha[i], dt.M[i])
        dt.Cx[i] = af.findCx(dt.Alpha[i], dt.M[i])
        dt.dy[i] = 0.5*dt.Cy[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dx[i] = 0.5*dt.Cx[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dt[i] = (dt.dy[i]*cos(dt.f[i])-dt.dx[i]*sin(dt.f[i]))*cos(dt.beta)
        dt.dq[i] = dt.dx[i]*cos(dt.f[i])+dt.dy[i]*sin(dt.f[i])
        dt.dmki[i] = dt.dy[i]*sin(dt.f[i])*r[i]*Iv/khi
        dt.dmkp[i] = dt.dx[i]*cos(dt.f[i])*r[i]*Cvt
        dt.G1[i] = Sigma*dt.Cy[i]*dt.W1[i]/8
        dST.append(dt.dt[i]*(r[i]-lhh))
        i = i + 1

    i = 1
    while i < N:
        dt.dY[i] = (dt.dy[i]+dt.dy[i-1])*dr[i-1]/2
        dt.ST = dt.ST + (dST[i]+dST[i-1])*dr[i-1]/2
        dt.dX[i] = (dt.dx[i]+dt.dx[i-1])*dr[i-1]/2
        dt.dT[i] = (dt.dt[i]+dt.dt[i-1])*dr[i-1]/2
        dt.dQ[i] = (dt.dq[i]+dt.dq[i-1])*dr[i-1]/2
        dt.dMki[i] = (dt.dmki[i]+dt.dmki[i-1])*dr[i-1]/2
        dt.dMkp[i] = (dt.dmkp[i] + dt.dmkp[i - 1]) * dr[i - 1] / 2
        dt.T = dt.T + dt.dT[i]
        dt.Mki = dt.Mki + dt.dMki[i]
        dt.Mkp = dt.Mkp + dt.dMkp[i]
        i = i + 1
    A1 = (Ihh*cos(dt.beta)-lhh*Shh)*(wR/R)**2*sin(dt.beta)
    dt.d2beta = (dt.ST - A1 - 9.81 * Shh) / (Ihh * (wR / R) ** 2)
    #print(dt.d2beta)
    print(degrees(dt.beta))
    dt.T = dt.T * z * khi
    dt.Mki = dt.Mki * z
    dt.Mkp = dt.Mkp * z
    dt.Mk = dt.Mki + dt.Mkp
    return dt


def twists(twist):
        k = []
        v = []
        for key, val in twist.items():
            k.append(key)
            v.append(val)
        i = 0
        j = 0
        while i < len(r):
            if k[j] <= r[i]:
                j = j + 1
            if j >= len(k):
                j = len(k)-1
            x = r[i]
            x1 = k[j-1]
            x2 = k[j]
            y1 = v[j-1]
            y2 = v[j]
            dFi.append((x-x1)*(y2-y1)/(x2-x1) + y1)
            i = i + 1


def th_calc(Fi7, dt, lN, hN, psi = 0):
    dt.Vyt = dt.Vy / wR
    dt.Fi7 = Fi7
    i = int(lN)
    rho = rho0*(20000-dt.H)/(20000+dt.H)
    dt.V_y = dt.Vy/wR
    while i < hN:
        dt.G0[i] = (dt.G0[i]+dt.G1[i])/2
        dt.u1[i] = dt.G0[i]/_r_[i]
        dt.v1[i] = -dt.V_y/2 + sqrt((dt.V_y/2)**2+abs(dt.G0[i]))*dt.G0[i]/abs(dt.G0[i])
        dt.V1[i] = dt.v1[i] + dt.V_y
        dt.U1[i] = _r_[i] - dt.u1[i]
        dt.W1[i] = sqrt(dt.V1[i]**2+dt.U1[i]**2)
        dt.f[i] = atan(dt.V1[i]/dt.U1[i])
        dt.fd[i] = degrees(dt.f[i])
        dt.Fi[i] = dt.Fi7+dFi[i]-K*dt.beta
        dt.Alpha[i] = dt.Fi[i] - dt.fd[i]
        dt.M[i] = dt.W1[i]*wR/340
        dt.Cy[i] = af.findCy(dt.Alpha[i], dt.M[i])
        dt.Cx[i] = af.findCx(dt.Alpha[i], dt.M[i])
        dt.dy[i] = 0.5*dt.Cy[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dx[i] = 0.5*dt.Cx[i]*rho*(dt.W1[i]*wR)**2*dt.b[i]
        dt.dt[i] = (dt.dy[i]*cos(dt.f[i])-dt.dx[i]*sin(dt.f[i]))*cos(dt.beta)
        dt.dq[i] = dt.dx[i]*cos(dt.f[i])+dt.dy[i]*sin(dt.f[i])
        dt.dmki[i] = dt.dy[i]*sin(dt.f[i])*r[i]*Iv/khi
        dt.dmkp[i] = dt.dx[i]*cos(dt.f[i])*r[i]*Cvt
        dt.G1[i] = Sigma*dt.Cy[i]*dt.W1[i]/8

        i = i + 1
    print("1--- %s seconds ---" % (time.time() - start_time))
    return dt


def calc_forces(dt):
    i = 1
    while i < N:
        dt.dY[i] = (dt.dy[i]+dt.dy[i-1])*dr[i-1]/2
        dt.dX[i] = (dt.dx[i]+dt.dx[i-1])*dr[i-1]/2
        dt.dT[i] = (dt.dt[i]+dt.dt[i-1])*dr[i-1]/2
        dt.dQ[i] = (dt.dq[i]+dt.dq[i-1])*dr[i-1]/2
        dt.dMki[i] = (dt.dmki[i]+dt.dmki[i-1])*dr[i-1]/2
        dt.dMkp[i] = (dt.dmkp[i] + dt.dmkp[i - 1]) * dr[i - 1] / 2
        dt.T = dt.T + dt.dT[i]
        dt.Mki = dt.Mki + dt.dMki[i]
        dt.Mkp = dt.Mkp + dt.dMkp[i]
        i = i + 1
    i = 0
    dt.T = dt.T * z * khi
    dt.Mki = dt.Mki * z
    dt.Mkp = dt.Mkp * z
    dt.Mk = dt.Mki + dt.Mkp


def blade_calc(blade_mass):
    global Shh
    global Ihh
    k = []
    v = []
    js = json.loads(blade_mass)
    rt = dict()
    for q in js[1:]:
        rt.update({q[0] : q[1]})
    for key, val in rt.items():
        k.append(key)
        v.append(val)
    Nm = len(k)
    Shh = 0
    Ihh = 0
    for i in range(Nm-1):
        Shh = Shh + (v[i] + v[i + 1]) / 2 * (k[i + 1] - k[i]) * (k[i] + k[i + 1]) / 2
        Ihh = Ihh + (v[i] + v[i + 1]) / 2 * (k[i + 1] - k[i]) * ((k[i] + k[i + 1]) / 2)**2


if __name__ == "__main__":
    print(degrees(atan2(0,0)))
    main()
