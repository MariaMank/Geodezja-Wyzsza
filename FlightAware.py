import math
import numpy as np
import csv
import plotly.graph_objects as go

a = 6378137 #metry
e2 = 0.00669438002290

def geo2xyz(fi, lam, h, a, e2):
    #fi = np.deg2rad(fi)
    #lam = np.deg2rad(lam)
    N = a / (math.sqrt(1 - e2*np.sin(fi)**2))
    x = (N+h)*np.cos(fi)*np.cos(lam)
    y = (N+h)*np.cos(fi)*np.sin(lam)
    z = (N*(1-e2)+h)*np.sin(fi)
    return x, y, z

def geo2neu(F1, L1, H1, F2, L2, H2 ):
    F1 = np.deg2rad(F1)
    L1 = np.deg2rad(L1)
    F2 = np.deg2rad(F2)
    L2 = np.deg2rad(L2)
    R = np.array([[-np.sin(F1)*np.cos(L1), -np.sin(L1), np.cos(F1)*np.cos(L1)],
                  [-np.sin(F1)*np.sin(L1), np.cos(L1), np.cos(F1)*np.sin(L1)],
                  [np.cos(F1), 0, np.sin(F1)]])
    x1, y1, z1 = geo2xyz(F1, L1, H1, a, e2)
    x2, y2, z2 = geo2xyz(F2, L2, H2, a, e2)
    x = np.array([[x2 - x1, y2 - y1, z2 - z1]])
    neu = np.dot(R.transpose(), x.transpose())
    return neu

#airport coordinates
f1 =35.341846
l1 = 25.148254
h1 = 35

lt = [f1]
lg = [l1]
nt = [0]
ea = [0]
up = [0]
result = open('NEUresult.txt', "w")
result.write("Date\t\tn\t\t\te\t\tu\t\tazymut\t\todleglosc\tkat z\n")
with open('FlightAware.txt') as data:
    file = csv.reader(data, delimiter='\t')
    line_count = 0
    czyzniknal = False
    for row in file:
        if line_count != 0:
            f2 = float(row[1])
            l2 = float(row[2])
            h2 = float(row[3])
            lt.insert(line_count, f2)
            lg.insert(line_count, l2)
            neu = geo2neu(f1, l1, h1, f2, l2, h2)
            n = neu[0]
            e = neu[1]
            u = neu[2]
            nt.insert(line_count, *n)
            ea.insert(line_count, *e)
            up.insert(line_count, *u)
            A = np.rad2deg(np.arctan(e / n))
            s = math.sqrt(n ** 2 + e ** 2 + u ** 2)
            z = np.rad2deg(np.arccos(u / s))
            line = row[0]+'\t'+str(*(np.round(n,3)))+'\t'+str(*(np.round(e,3)))+'\t'+str(*(np.round(u,3)))+'\t'+\
                   str(*(np.round(A, 6)))+'\t'+str((np.round(s,3)))+'\t'+str(*(np.round(z,6)))+'\n'
            result.write(line)
            if czyzniknal == False and (z >= 90 and u <= 0):
                print('Samolot zniknął za horyzontem w punkcie o współrzędnych neu: [', *(np.round(n,3)), *(np.round(e,3)), *(np.round(u,3)), '] i kącie zenitalnym: ',
                      *(np.round(z,6)),'\nI współrzędnych geodezyjnych fi, lamda, h: [', row[1], row[2], row[3], ']')
                fig = go.Figure(go.Scattermapbox(
                    mode="markers",
                    lon=[l2],
                    lat=[f2],
                    name= 'Miejsce gdzie samolot znika za horyzontem',
                    marker=go.scattermapbox.Marker(
                        size=6,
                        color='rgb(0, 0, 0)'
                    )))
                czyzniknal = True
            line_count += 1
        else:
            line_count += 1
result.close()

if czyzniknal!=True:
    print("Samolot nie zniknął za horyzontem")
    fig = go.Figure(go.Scattermapbox(
        mode="lines",
        lon=lg,
        lat=lt,
        name='Trasa samolotu',
        marker=go.scattermapbox.Marker(
                size=8,
                color='rgb(192, 20, 20)'
        )))
else:
    fig.add_trace(go.Scattermapbox(
        mode="lines",
        lon=lg,
        lat=lt,
        name='Trasa samolotu',
        marker=go.scattermapbox.Marker(
                size=8,
                color='rgb(192, 20, 20)'
        )))

fig.update_layout(
    margin={'autoexpand': False, 'l': 0, 't': 45, 'b': 0, 'r': 0},
    title={
        'font': {'color': 'rgb(0, 0, 0)', 'family': "Arial", 'size':30},
        'x': 0.5,
        'y': 0.987,
        'text': 'Trasa lotu samolotu TUI3GB Heraklion - Hanover 28.10.2021'
        },
    hovermode='closest',
    mapbox={
        'style': "stamen-terrain",
        'center': {'lon': 18, 'lat': 44.554},
        'zoom': 4.555
       },
    legend={
        'title': 'Legenda',
        'font': {'family': "Arial"},
        'bordercolor': "grey",
        'borderwidth': 5,
        'y': 0,
        'x': 0.85
        },
    paper_bgcolor="lightgrey"
)


fig3d = go.Figure(data=[go.Scatter3d(
    x=nt,
    y=ea,
    z=up,
    mode='markers',
    marker=dict(
        size=8,
        color="red",
        opacity=0.8
    ),
)])

fig3d.update_layout(margin=dict(l=0, r=0, b=0, t=0),
                    scene={
                        "xaxis": {"title": 'n'},
                        "yaxis": {"title": 'e'},
                        "zaxis": {"title": 'u'},
                        'camera_eye': {"x": -1, "y": -2, "z": 1}
                    },
                    title={
                        'font': {'color': 'rgb(0, 0, 0)', 'family': "Arial", 'size': 30},
                        'x': 0.5,
                        'y': 0.987,
                        'text': 'Trasa lotu samolotu TUI3GB Heraklion - Hanover 28.10.2021 we współrzednych neu'
                    },)
fig3d.show()
fig.show()


print("Przeliczone współrzędne neu wraz z Azymutami, odlełością skośna i kątem zenitalnym zostały zapisane do pliku NEUresult.txt")
