import os
from main import aprox

Alpha = []
Cy = []
Cx = []


def read():
    dir = ".\\NACA 63012\\"
    flist = os.listdir(dir)
    ofile = []
    uM = dict()
    global Alpha
    global Cy
    global Cx
    global M
    i = 0
    for file in flist:
        fname = dir+file
        ofile.append(open(fname, 'r'))
        rfile = ofile[i].readlines()
        uM.update({fname : float(rfile[7][10:15])})
        i = i + 1
    M = sorted(uM.items(), key=lambda x : x[1])
    i = 0
    for Q in M:
        Alpha.append([])
        Cy.append([])
        Cx.append([])
        fname = Q[0]
        ofile1 = open(fname, 'r')
        rfile = ofile1.readlines()
        for L in rfile[11:-2]:
            Alpha[i].append(float(L[0:7]))
            Cy[i].append(float(L[10:17]))
            Cx[i].append(float(L[20:27]))
        i = i + 1


def findCy(A, m = 0.1):
    i = 0
    while i < len(M):
        if m <= M[i][1]:
            break
        i = i + 1
    x1 = M[i-1][1]
    x2 = M[i][1]
    y1 = aprox(A, Alpha[i-1], Cy[i-1])
    y2 = aprox(A, Alpha[i], Cy[i])
    return (m - x1) * (y2 - y1) / (x2 - x1) + y1


def findCx(A, m = 0.1):
    i = 0
    while i < len(M):
        if m <= M[i][1]:
            break
        i = i + 1
    x1 = M[i-1][1]
    x2 = M[i][1]
    y1 = aprox(A, Alpha[i-1], Cx[i-1])
    y2 = aprox(A, Alpha[i], Cx[i])
    return (m - x1) * (y2 - y1) / (x2 - x1) + y1
