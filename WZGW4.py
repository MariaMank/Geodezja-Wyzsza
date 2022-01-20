import numpy as np
import shapely.geometry as shp

a = 6378137
e2 = 0.00669438002290
punktyfi = [50.25, 50.0, 50.25 ,50.0 ,50.125, 50.125270449027546]        #punkty wejściowe
punktylm = [20.75, 20.75, 21.25, 21.25 ,21.0, 21.00065108883011]
p2000 = []                      #tablice wyników xy, xy itd
p1992 = []
pgk2 = []
pgk0 = []
mk = []
m2k2 = []
def Np(f, a, e2):
    N = a / np.sqrt(1 - e2 * (np.sin(f) ** 2))
    return N
def Mp(f, a, e2):
    M = (a * (1 - e2)) / ((np.sqrt(1 - e2 * np.sin(f) ** 2)) ** 3)
    return M
def m(y, r):
    mp = 1 + (y**2)/(2*(r**2)) + (y**4)/(24*(r**4))
    return mp

def geo2GK(fi, lm, lm0): #coś nie tak
    fi = np.deg2rad(fi)
    lm = np.deg2rad(lm)
    lm0 = np.deg2rad(lm0)
    ep2 = 0.0067394909677548 #kwadrat drugiego mimośrodu
    e2 = 0.006694380022900  #kwarat pierwszego mimośrodu
    a = 6378137
    N = a/np.sqrt(1 - e2*(np.sin(fi))**2)
    t = np.tan(fi)
    n2 = ep2*(np.cos(fi)**2)
    l = lm - lm0
    a0 = 1 - (e2/4) - (3*e2**2)/64 - (5*e2**3)/256
    a2 = (3/8)* ((e2 + (e2**2)/ 4) + (15 * (e2 ** 3)) / 128)
    a4 = (15/256)*(e2**2+(3*e2**3)/4)
    a6 = (35*(e2**3))/3072
    sigma = a*(a0*fi - a2*np.sin(2*fi)+ a4*np.sin(4*fi) - a6*np.sin(6*fi))
    x = sigma + ((l**2)/2)*N*np.sin(fi)*np.cos(fi)*(1+((l**2)/12)*(np.cos(fi)**2)*(5 - t**2 + 9*n2 + 4*(n2)**2) \
        + ((l**4)/360)*(np.cos(fi)**4)*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*(t**2)))
    y = l*N*np.cos(fi)*(1+((l**2)/6)*(np.cos(fi)**2)*(1 - t**2+n2) +
                        ((l**4)/120)*(np.cos(fi)**4)*(5 - 18*(t**2) + t**4 + 14*n2 - 58*n2*(t**2)))
    return x, y

def gk2u92 (xk, yk):
    m = 0.9993
    x = xk*m - 5300000
    y = yk*m + 500000
    return x, y
def l(lm):
    if lm >22.5 and lm < 25.5:
        lm0 = 24
    elif lm >19.5 and lm < 22.5:
        lm0 = 21
    elif lm < 19.5 and lm > 16.5:
        lm0 = 18
    elif lm < 13.5 and lm < 16.5:
        lm0 = 15
    return lm0
def gk2u00 (xk, yk, lm):
    m = 0.999923
    lm0 = l(lm)
    c = lm0/3
    x = m*xk
    y = m*yk + 500000 + c*1000000
    return x, y

def u922gk (x, y):
    m = 0.9993
    xk = (x + 5300000)/m
    yk = (y - 500000)/m
    return xk, yk


def u002gk (x, y,lm0):
    m = 0.999923
    c = lm0/3
    xk = x/m
    yk = (y - c*1000000 - 500000)/m
    return xk, yk

def gk2el (xk, yk,lm0, a, e2):
    f = 1/(298.257)
    b = (a - (f*a))
    e22 = ((a**2 - b**2)/b**2)
    ep = 1
    lm0 = np.deg2rad(lm0)
    a0 = 1 - (e2 / 4) - (3 * e2 ** 2) / 64 - (5 * e2 ** 3) / 256
    a2 = (3 / 8) * ((e2 + (e2 ** 2) / 4) + (15 * (e2 ** 3)) / 128)
    a4 = (15 / 256) * (e2 ** 2 + (3 * e2 ** 3) / 4)
    a6 = (35 * (e2 ** 3)) / 3072
    fi0 = xk/(a*a0)
    sigma = a*(a0*fi0 - a2*np.sin(2*fi0)+ a4*np.sin(4*fi0) - a6*np.sin(6*fi0))
    eps = np.deg2rad(0.00001/3600)
    while ep> eps:
        fi1 = fi0 + (xk - sigma)/(a* a0)
        if abs(fi1 - fi0) < eps:
            break
        else:
            fi0 = fi1
        sigma = a * (a0 * fi1 - a2 * np.sin(2 * fi1) + a4 * np.sin(4 * fi1) - a6 * np.sin(6 * fi1))

    N = Np(fi1, a, e2)
    M = Mp(fi1, a, e2)
    t = np.tan(fi1)
    n2 = e22*(np.cos(fi1)**2)
    fi = fi1 - (((yk**2)*t)/(2*M*N) )* (1 - (yk**2)/(12*(N**2))*(5 + 3*(t**2) + n2 - 9*n2*(t**2) - 4*(n2**2))
                                      + (yk**4)/(360*(N**4))*(61 + 90 * (t**2) + 45*(t**4)))
    lm = lm0 + (yk/(N*np.cos(fi1))) * (1 - (yk**2)/(6*(N**2)) * (1 + 2 * (t**2)+n2)
         + ((yk**4)/(120*(N**4)))* (5 + 28* (t**2) + 24*(t**4) + 6* n2 + 8*n2*(t**2)))
    fi = np.rad2deg(fi)
    lm = np.rad2deg(lm)
    return fi, lm

def mkappa (x, y, a, e2, u, lmm, fi):
    if u == '92':
        m0 = 0.9993
        gk = u922gk(x, y)
        yk = gk[1]
        N = Np(fi, a, e2)
        M = Mp(fi, a, e2)
        R = np.sqrt(M*N)
        m1 =m(yk, R)*m0
        m2 = (m(yk, R)**2)*(m0**2)
        k = (m1 - 1)*1000
        k2 = (m2 - 1)*10000
    elif u == '00':
        m0 = 0.999923
        gk = u002gk(x, y, 21)
        yk = gk[1]
        N = Np(fi, a, e2)
        M = Mp(fi, a, e2)
        R = np.sqrt(M*N)
        m1 =m(yk, R)*m0
        m2 = (m(yk, R)**2)*(m0**2)
        k = (m1 - 1)*1000
        k2 = (m2 - 1)*10000
    elif u == 'None':
        N = Np(fi, a, e2)
        M = Mp(fi, a, e2)
        R = np.sqrt(M*N)
        m1 = m(y, R)
        m2 = m1**2
        k = (m1 - 1)*1000
        k2 = (m2 - 1) * 10000
    return m1, k, m2, k2

def pole(a, e2, lm1, lm2, fi1, fi2): #pole ze współrzędnych geo z zad 3
    fi1 = np.deg2rad(fi1)
    lm1 = np.deg2rad(lm1)
    fi2 = np.deg2rad(fi2)
    lm2 = np.deg2rad(lm2)
    b = a * np.sqrt(1 - e2)
    e = np.sqrt(e2)
    a1 = (np.sin(fi1)/(1 - e2*np.sin(fi1)**2) + (1/(2*e))*np.log((1 + e*np.sin(fi1))/(1 - e*np.sin(fi1))))
    a2 = (np.sin(fi2)/(1 - e2*np.sin(fi2)**2) + (1/(2*e))*np.log((1 + e*np.sin(fi2))/(1 - e*np.sin(fi2))))
    p = b**2*(lm2 - lm1)/2*(a1 - a2)
    return p

def wywolywanie():
    p = 0
    print('------------------------------ZESTAWIENIE WSPÓŁRZĘDNYCH---------------------------------')
    line = 'pkt' + '\t' + 'fi' + '\t\t' + 'lambda' + '\t' + 'X GK00\t\t' + 'Y GK00' + '\t\t' + 'X GK92' + \
           '\t\t' + 'Y GK92' + '\t\t' + 'X 2000' + '\t\t' + 'Y 2000' + '\t\t' + 'X 1992' + '\t\t' + 'Y 1992'
    print(line)
    for i in punktyfi:
        fi = i
        lm = punktylm[p]
        lm0 = l(lm)
        gkk = geo2GK(fi, lm, lm0)       #GK dla 00  - odpowiedni poludnik osiowy
        gkkk = geo2GK(fi, lm, 19)       #GK dla 92 - poludnik osiowy 19
        pgk0.append([gkk[0], gkk[1]])
        pgk2.append([gkkk[0], gkkk[1]])
        u200 = gk2u00(gkk[0], gkk[1], lm0) #GK00 ---> 2000
        m00 = mkappa(u200[0], u200[1], a, e2, '00', lm0, fi)
        p2000.append([u200[0], u200[1]])
        u1992 = gk2u92(gkkk[0], gkkk[1])   #GK92---> 1992
        m92 = mkappa(u1992[0], u1992[1], a, e2, '92', 19, fi)
        p1992.append([u1992[0], u1992[1]])
        u92kk = u922gk(u1992[0], u1992[1]) #92 --->GK92
        u20kk = u002gk(u200[0], u200[1], lm0) #2000 --->GK00
        mGK0 = mkappa(u20kk[0], u20kk[1], a, e2, 'None', lm0, fi)
        mGK2 = mkappa(u92kk[0], u92kk[1], a, e2, 'None', 19, fi)
        gkgeo = gk2el(u92kk[0], u92kk[1], 19, a, e2) #GK92---> fi lm GRS80
        gkgeo1 = gk2el(u20kk[0], u20kk[1], lm0, a, e2) #GK00 ---> fi lm GRS80
        mk.append([mGK0[0],mGK0[1], mGK2[0], mGK2[1], m00[0], m00[1], m92[0], m92[1]])
        m2k2.append([mGK0[2], mGK0[3], mGK2[2], mGK2[3], m00[2], m00[3], m92[2], m92[3]])
        line = str(p + 1) + '\t' + str(fi) + '\t' + str(lm) + '\t' + str(round(gkk[0], 3))+'\t'+str(round(gkk[1], 3))+ '\t' \
               + str(round(gkkk[0], 3))+ '\t'+str(round(gkkk[1], 3)) + '\t' + str(round(u200[0], 3))+ '\t'+str(round(u200[1], 3)) + \
               '\t'+ str(round(u1992[0], 3))+ '\t'+str(round(u1992[1], 3))
        print(line)
        p+=1
    print('------------------------------ZESTAWIENIE PÓL POWIERZCHNI---------------------------------')
    pol = pole(a, e2, punktylm[0], punktylm[3], punktyfi[0], punktyfi[3]) / 1000000
    pp92 = shp.Polygon([(p1992[0][0], p1992[0][1]),(p1992[1][0], p1992[1][1]), (p1992[3][0], p1992[3][1]), (p1992[2][0], p1992[2][1])] )
    pol92 = pp92.area/1000000
    pol00 = (np.abs(p2000[3][0]  - p2000[0][0])*np.abs(p2000[0][1]  - p2000[3][1]))/ 1000000
    ppgko = shp.Polygon([(pgk0[0][0], pgk0[0][1]), (pgk0[1][0], pgk0[1][1]), (pgk0[3][0], pgk0[3][1]),
                         (pgk0[2][0],pgk0[2][1])])
    polgk = ppgko.area / 1000000
    ppgk2 = shp.Polygon([(pgk2[0][0], pgk2[0][1]), (pgk2[1][0], pgk2[1][1]), (pgk2[3][0], pgk2[3][1]),
                         (pgk2[2][0],pgk2[2][1])])
    polg2 = ppgk2.area / 1000000
    line1 = 'pole elip [km^2]' + '\t' + 'pole GK92 [km^2]' + '\t' + 'pole GK00 [km^2]\t' + 'pole 2000 [km^2]' + '\t' + 'pole 1992 [km^2]'
    print(line1)
    line2 = str(round(pol, 12)) + '\t' +str(round(polg2, 12))   + '\t' + str(round(polgk, 12)) + '\t' + str(round(pol00, 12))  +'\t' + str(round(pol92, 12))
    print(line2)
    print('--------------------------ELEMENTARNA SKALA DŁUGOŚCI I ZNIEKSZTAŁCENIA NA KM---------------------------')
    line = 'pkt' + '\t' + 'mGK00' + '\t\t' + 'KGK00' + '\t' + 'mGK92\t\t' + 'KGK92' + '\t\t' + 'm00' + '\t\t' + 'K00' + '\t\t' + 'm92' + '\t\t' + 'K92'
    print(line)
    ke = 0
    for i in punktyfi:
        print(ke, *mk[ke])
        ke+=1
    print('--------------------------ELEMENTARNA SKALA DŁUGOŚCI I ZNIEKSZTAŁCENIA NA HA---------------------------')
    line = 'pkt' + '\t' + 'mGK00^2' + '\t\t' + 'KGK00^2' + '\t' + 'mGK92^2\t\t' + 'KGK92^2' + '\t\t' + 'm00^2' + '\t\t' + 'K00^2' + '\t\t' + 'm92^2' + '\t\t' + 'K92^2'
    print(line)
    ke = 0
    for i in punktyfi:
        print(ke, *m2k2[ke])
        ke+=1

wywolywanie()