import numpy as np
#----PARAMETRY GRS80
a = 6378137 #metry
e2 = 0.00669438002290
#---PARAMETRY KRASOWSKI
a2 = 6378245
e22 = (0.0818133340169)**2

punktyfi = [50.25, 50.0, 50.25 ,50.0 ,50.125, 50.125270449027546]        #punkty wejściowe
punktylm = [20.75, 20.75, 21.25, 21.25 ,21.0, 21.00065108883011]

def xyzg2xyzk (xg, yg, zg):
    Tx = -33.4297
    Ty = 146.5746
    Tz = 76.2865
    m = 1 + 0.8407728 /10**6
    ex = -1.7388854 / 10**6
    ey = -0.2561460 * 10**(-6)
    ez = 4.0896031 * 10**(-6)
    xk = m * xg + ez * yg - ey* zg + Tx
    yk = (-ez) * xg + m * yg + ex* zg + Ty
    zk = ey * xg - ex * yg + m* zg + Tz
    return xk, yk, zk

def rounding(x, a):
    p = 0
    for i in x:
        x[p] = round(i, a)
        p+=1
    return x

def Np(a, e2, fi):
    N = a / (np.sqrt(1 - e2*np.sin(fi)**2))
    return N

def Hirvonen(x, y, z, a, e2):
    ep = np.deg2rad(0.000001/3600)
    r = np.sqrt(x**2+y**2)
    f1= np.arctan(z/(r*(1-e2)))
    f2 = 2*f1
    while np.abs(f1 - f2) > ep:
        f2 = f1
        N =  Np(a, e2, f1)
        h = r / np.cos(f1) - N
        f1 = np.arctan(z/(r*(1-e2*N/(N+h))))
    lam = np.arctan(y/x)
    N1 = Np(a, e2, f1)
    h1 = r / np.cos(f1) - N
    f2 = np.rad2deg(f2)
    lam = np.rad2deg(lam)
    return f2, lam, N1, h1

def geo2xyz(fi, lam ,h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = Np(a, e2, fi)
    x = (N+h)*np.cos(fi)*np.cos(lam)
    y = (N+h)*np.cos(fi)*np.sin(lam)
    z = (N*(1-e2)+h)*np.sin(fi)
    return x, y, z

def flg2flk():
    p = 0
    line = 'pkt' +'\t'+ 'fi G..' + '\t' + 'lam G' + '\t' + 'H G\t' + 'X GRS80' + '\t\t' + 'Y GRS80' + '\t\t' + 'Z GRS80'+ '\t\t' + 'X Krasowski' + '\t' + 'Y Krasowski' + '\t' + 'Z Krasowski' + '\t' + 'fi Kras...' + '\t' + 'lambda Kras' + '\t' + 'H Kras'
    print(line)
    for i in punktyfi:
        fi = i
        lm = punktylm[p]
        N = Np(a, e2, fi)
        h = 100
        xg, yg, zg = geo2xyz(fi, lm, h, a, e2)            #fi, lamba GRS80 na xyz na GRS
        xk, yk, zk = xyzg2xyzk(xg, yg, zg)                   #xyz na GRS na xyz na Krasowkim
        fik, lmk, Nk, hk = Hirvonen(xk, yk, zk, a2, e22)  #xyz Krasowski na fi lambda krasowski
        xg, yg, zg = rounding([xg, yg, zg], 3)  # fi, lamba GRS80 na xyz na GRS
        xk, yk, zk = rounding([xk, yk, zk ], 3)  # xyz na GRS na xyz na Krasowkim
        fik, lmk  = rounding([fik, lmk], 6)
        hk = rounding([hk], 3)
        line = str(p+1) + '\t'+ str(fi) +'\t'+ str(lm) +'\t'+  str(h) + '\t'+str(xg)+'\t'+  str(yg) +'\t'+ str(zg) +'\t'+  str(xk) + '\t'+ str(yk) +'\t'+  str(zk) + '\t'+ str(fik) +'\t'+  str(lmk) +'\t'+  str(*hk)
        p+=1
        print(line)


flg2flk()