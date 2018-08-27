class Data:

    def __init__(self, n, Fi7, psi = 0, Vy = 0, Vx = 0, Vz = 0):
        i = 0
        self.Fi7 = Fi7
        self.H = 0
        self.AlphaNV = 0
        self.psi = psi
        self.dPsi = 0
        self.G0 = []
        self.u1 = []
        self.Vy = Vy
        self.Vx = Vx
        self.V_y = 0
        self.v1 = []
        self.V1 = []
        self.U1 = []
        self.W1 = []
        self.beta = 0
        self.f = []
        self.fd = []
        self.Fi = []
        self.Alpha = []
        self.M = []
        self.Cy = []
        self.Cx = []
        self.b = []
        self.dy = []
        self.dY = []
        self.dx = []
        self.dX = []
        self.dt = []
        self.dT = []
        self.T = 0
        self.Mki = 0
        self.Mkp = 0
        self.Mk = 0
        self.dq = []
        self.dQ = []
        self.dmki = []
        self.dMki = []
        self.dmkp = []
        self.dMkp = []
        self.G1 = []
        self.ST = 0
        self.d2beta = 0
        self.dbeta = 0
        while i < n:
            self.M.append(0)
            self.G0.append(0.004)
            self.u1.append(0)
            self.v1.append(0)
            self.V1.append(0)
            self.U1.append(0)
            self.W1.append(0)
            self.f.append(0)
            self.fd.append(0)
            self.Fi.append(0)
            self.Alpha.append(0)
            self.Cy.append(0)
            self.Cx.append(0)
            self.b.append(0.17)   #Заглушка!
            self.dy.append(0)
            self.dY.append(0)
            self.dx.append(0)
            self.dX.append(0)
            self.dt.append(0)
            self.dT.append(0)
            self.dq.append(0)
            self.dQ.append(0)
            self.dmki.append(0)
            self.dMki.append(0)
            self.dmkp.append(0)
            self.dMkp.append(0)
            self.G1.append(0)
            i = i + 1

    def add_r(self, i, it):
        self.r[i] = it
        print(self.r[i])


