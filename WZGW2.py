import numpy as np
import jdcal as jc
import tkinter as tk, tkinter.ttk as ttk
import matplotlib.backends.backend_tkagg as TkAgg
from matplotlib.figure import Figure

window = tk.Tk()
window.geometry("900x800")
window.title("ETP - projekt1")
window.config(background="#c0d6e4")
tabControl = ttk.Notebook(window)
tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)
tab3 = ttk.Frame(tabControl)
tabControl.add(tab1, text='półkula północna')
tabControl.add(tab2, text='okolice równika')
tabControl.add(tab3, text='półkula południowa')
tabControl.pack(expand=1, fill="both")

#Denebola - lew
a = 11.81777778 #rektascensja
#a = 168.2711111 #alfa
dek = np.deg2rad(14.57672222) #deklinacja
fi1= np.deg2rad(45)
lb1 = np.deg2rad(17)
fi2= np.deg2rad(1)
lb2= np.deg2rad(25)
fi3= np.deg2rad(-60)
lb3= np.deg2rad(190)

#2. Przeliczenie czasu słonecznego UT na czas gwiazdowy S oraz obliczenie kąta godzinnego
def GMST (jd): #dobrze
    t = (jd - 2451545)/36525
    g = 280.46061837 + 360.98564736629*(jd - 2451545) + 0.000387933*(t**2) - (t**3)/38710000
    g = g%360
    return g

def kh (y, m, d, h, lb, alfa): #dobrze
    jd = sum(jc.gcal2jd(y, m, d))
    gm = GMST(jd)
    UT1 = h * 1.002737909350795
    S = UT1*15 + lb + gm #czas gwiazdowy (stopnie)
    t = S - alfa*15 # kat godzinny (stopnie)
    t = np.deg2rad(t) #kat godzinny radiany
    return t

def figury(h, A, z, tab, x, y, zz):
    fig = Figure(figsize=(5, 4), dpi=100)
    fig.add_subplot(111).plot(h, A)
    fig.supxlabel('godzina', size=10)
    fig.supylabel('Azymut [stopnie]', size=10)
    fig.suptitle('Zmiana azymutu w czasie')
    fig2 = Figure(figsize=(5, 4), dpi=100)
    fig2.add_subplot(111).plot(h, z)
    fig2.supxlabel('godzina', size=10)
    fig2.supylabel('odległość zenitalna', size=10)
    fig2.suptitle('Zmiana odelgłości zenitalnej w czasie')
    fig3 = Figure( figsize=(5,4),dpi = 100)

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    sf1 = np.outer(np.cos(u), np.sin(v))
    sf2 = np.outer(np.sin(u), np.sin(v))
    sf3 = np.outer(np.ones(np.size(u)), np.cos(v))
    ax = fig3.add_subplot(111, projection='3d')
    ax.plot_surface(sf1, sf2, sf3, alpha=0.3)
    ax.scatter3D(x, y, zz, c='yellow',alpha=1)
    fig3.suptitle('Położenie gwiazdy')
    canvas = TkAgg.FigureCanvasTkAgg(fig, master=tab)  # A tk.DrawingArea.
    canvas2 = TkAgg.FigureCanvasTkAgg(fig2, master=tab)
    canvas3 = TkAgg.FigureCanvasTkAgg(fig3, master=tab)
    canvas.draw()
    canvas2.draw()
    canvas3.draw()
    canvas3.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas2.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)

def A(fi, dek, t): #w poprawnych jednostkach już
    g = -np.cos(dek)*np.sin(t)
    d = np.cos(fi)*np.sin(dek)-np.sin(fi)*np.cos(dek)*np.cos(t)
    a = np.arctan(g/d)
    ad = np.rad2deg(a)
    if g>0:
        if d<0:
            ad = 180 + ad
        if d> 0 and ad<360:
            return a
    elif g<0:
        if d>0:
            ad= 360 + ad
        elif d<0:
            ad+=180
    if ad > 360:
        ad-=360
    elif ad < 0:
        ad+=360
    a = np.deg2rad(ad)
    return a #też w radianach

def z(fi, dek, t):
    res = np.arccos(np.sin(fi)*np.sin(dek)-np.cos(fi)*np.cos(dek)*np.cos(t))
    return res


def toxyz (r, z, A):
    x = r*np.sin(z)*np.cos(A)
    y = r*np.sin(z)*np.sin(A)
    z = r*np.cos(z)
    return x, y, z
def licz(fi, lb):
    x1 = []
    y1= []
    z1 = []
    od1 = []
    A1= []
    h = []
    for i in range(24):
        t = kh(2022, 8, 13, i, lb, a)
        d = z(fi, dek, t)
        Az = A(fi, dek, t)
        h.append(i)
        resl1 = toxyz(1, d, Az )
        od1.append(d)
        A1.append(np.rad2deg(Az))
        x1.append(resl1[0])
        y1.append(resl1[1])
        z1.append(resl1[2])
    return h, A1, od1, x1, y1, z1

def glowna():
    a = licz(fi1, lb1)
    a1= licz(fi2, lb2)
    a2 = licz(fi3, lb3)
    figury(a[0], a[1], a[2], tab1, a[3], a[4], a[5])
    figury(a1[0], a1[1], a1[2], tab2, a1[3], a1[4], a1[5])
    figury(a2[0], a2[1], a2[2], tab3, a2[3], a2[4], a2[5])

glowna()

window.mainloop()