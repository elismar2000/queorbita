# Programa "plot_queorbita.py
# Lê arquivo de saída do "queorbita"
# Pra carda órbita, roda o "pot_hbd_din_fric-2" pra frente e pra trás no tempo
#    por 1 milhão de anos. Os arquivos de saída ficam no diretório
# Plota cada órbita
#
#
# FALTA Botar argumentos de entrada:
#    - Quais órbitas plotar. Variável "plotar".
#        f=com dyn fric
#        nf=sem dyn fric
#        kepl=kepleriana
#        sintaxe atual:    plotar=["f","nf","kepl"]
#
#    - Seleção das órbitas que quer plotar
#        sintaxe atual:
#          a=queorb.loc[(queorb["e"] == 0.9) & (queorb["q"] == 5.7) & (queorb["spy"] >= 0.8) & (queorb["PERIC"] == "POS") & (queorb["Dir"] == "Pro")]
#          orb_N=a.index
#
#
import numpy as np
import pandas as pd                 # Para manipular arquivos: read_csv
import plotly.graph_objects as go   # Para plot em 3d
import plotly.express as px
import subprocess                   # Para rodar programas externos
import os                           # Para rodar programas externos
import sys


queorb = pd.read_csv('queorbita_new.out' ,  header=53, sep="\t")

#queorb.dtypes

#queorb.head()

#queorb.tail()

#print(queorb)

e           = queorb['e']
q           = queorb['q']
rnow        = queorb['rnow']
vq          = queorb['vq']
gx          = queorb['gx']
gy          = queorb['gy']
gz          = queorb['gz']
VX_ceu      = queorb['VX_ceu']
VY_ceu      = queorb['VY_ceu']
VZ_ceu      = queorb['VZ_ceu']
x           = queorb['x']
y           = queorb['y']
vx          = queorb['vx']
vy          = queorb['vy']
spx         = queorb['spx']
spy         = queorb['spy']
spz         = queorb['spz']
qx          = queorb['qx']
qy          = queorb['qy']
qz          = queorb['qz']
vsys        = queorb['vsys']
Dir         = queorb['Dir']
# spin-orb    = queorb['spin-orb']
PERIC       = queorb['PERIC']
# pos-peri    = queorb['pos-peri']



# Crio diretório para arquivos temporários

dirName = "orbits_temp"
if os.path.exists(dirName):
    print("Directory " , dirName ,  " exists.  ")
else:
    sys.exit("Directory " + dirName +  " does not exist. Please run calcula_orbitas_queorbita.py first.")

# Começo a povoar o gráfico com o marcador do centro do sistema de coordenadas:

# Marcador do centro do sistema de coordenadas
c = 0.0
fig = go.Figure(data=[go.Scatter3d(
    x=[c], y=[c], z=[c],
    mode='markers',
    marker_size=2,
    name='centro',
#    showlegend=False,
#    opacity=0.2
    marker=dict(
        size=2,
        color='rgb(0, 0, 0)',        # set color to an array/list of desired values
#        opacity=0.8
    ),
)])


# Escolho as órbitas que serão plotadas:
##################### Orbita 1
#orb_N = np.arange(1000,1100,1)
orb_N = [0:9]    # No arquivo do queorbita --> numero da linha menos 53
                                      # aqui Orbita_07 Orbita_08 Orbita_09 Orbita_10 Orbita_11
#orb_N = [2799, 2801, 2803, 2805, 2807, 2809, 2811, 2813, 2814, 2816, 2822, 2823, 2824, 2825] # orbitas com e=1 q=2.5
a=queorb.loc[(queorb["e"] == 0.9) & (queorb["q"] == 5.7) & (queorb["spy"] >= 0.8) & (queorb["PERIC"] == "POS") & (queorb["Dir"] == "Pro")]
orb_N=a.index
orb_N = [2555, 5987]

# Qual orbita plotar?
#   f=com dyn fric
#   nf=sem dyn fric
#   kepl=kepleriana
#
plotar=["f","nf","kepl"]
#plotar=[nf"]
#plotar=[kepl"]
#plotar=["f","nf"]
#plotar=["f"]


# Loop sobre todas as órbitas pra desenhar cada uma delas:
for i in orb_N:
    fname = dirName + "/" + "orb_%05d.dat" % i
    orb = pd.read_csv(fname, sep="\t")
    orb1 = orb[9000:11000:10]       # Pega somente 1Gyr antes e 1Gyr depois do pericentro (linha 10000).
                                    # Os índices ficam os mesmos do dataframe antigo
    orb1 = orb1.reset_index(drop=True)     # Arruma os índices do orb1 em sequência a partir do 0

# Com dynamical friction
    if "f" in plotar:
        fig.add_trace(go.Scatter3d(
            x=orb1.x_f,
            y=orb1.y_f,
            z=orb1.z_f,
            name="orb_%05d fric" % i,
            mode='lines',
#    marker=dict(
#        size=2,
#        color=orb1.time,        # set color to an array/list of desired values
#        colorscale='bluered',   # choose a colorscale
#        opacity=0.8
#    ),
            line=dict(
                color=i,
#        color=orb1.time,
                width=2
            )
        )
                     )
# Marcador do início da órbita COM dynamical friction
        x_c, y_c, z_c = orb1.x_f[0], orb1.y_f[0], orb1.z_f[0]
        fig.add_trace(go.Scatter3d(
            x=[x_c],
            y=[y_c],
            z=[z_c],
            mode='markers',
            marker_size=2,
            name='início',
            showlegend=False,
#            marker_color=0,
#            opacity=0.8,
            marker=dict(
                size=2,
                color='rgb(0, 0, 255)',        # set color to an array/list of desired values
                opacity=0.8
            ),
        )
                     )
# FIM Com dynamical friction





# Sem dynamical friction
    if "nf" in plotar:
        fig.add_trace(go.Scatter3d(
            x=orb1.x_nf,
            y=orb1.y_nf,
            z=orb1.z_nf,
            name="orb_%05d" % i,
            mode='lines',
#    marker=dict(
#        size=2,
#        color=orb1.time,        # set color to an array/list of desired values
#        colorscale='bluered',   # choose a colorscale
#        opacity=0.8
#    ),
            line=dict(
                color=i,
#        color=orb1.time,
                width=2
            )
        )
                     )
# Marcador do início da órbita SEM dynamical friction
        x_c, y_c, z_c = orb1.x_nf[0], orb1.y_nf[0], orb1.z_nf[0]
        fig.add_trace(go.Scatter3d(
            x=[x_c],
            y=[y_c],
            z=[z_c],
            mode='markers',
            marker_size=2,
            name='início',
            showlegend=False,
            marker_color=0,
            opacity=0.8
        )
                     )
# FIM Sem dynamical friction



# Kepleriana
    if "kepl" in plotar:
        fig.add_trace(go.Scatter3d(
            x=orb1.x_kepl,
            y=orb1.y_kepl,
            z=orb1.z_kepl,
            name="orb_%05d kepler" % i,
            mode='lines',
#    marker=dict(
#        size=2,
#        color=orb1.time,        # set color to an array/list of desired values
#        colorscale='bluered',   # choose a colorscale
#        opacity=0.8
#    ),
            line=dict(
                color=i,
#        color=orb1.time,
                width=2
            )
        )
                     )

# Marcador do início da órbita Kepleriana
        x_c, y_c, z_c = orb1.x_kepl[0], orb1.y_kepl[0], orb1.z_kepl[0]
        fig.add_trace(go.Scatter3d(
            x=[x_c],
            y=[y_c],
            z=[z_c],
            mode='markers',
            marker_size=2,
            name='início',
            showlegend=False,
            marker_color=0,
            opacity=0.8
        )
                     )
# FIM Kepleriana







# Marcador do pericentro
    a=queorb.loc[queorb.index==i]
#    print(a.qx[i],a.qy[i],a.qz[i],a.q[i])
    x_p=[0.0,a.qx[i]*a.q[i]]
    y_p=[0.0,a.qy[i]*a.q[i]]
    z_p=[0.0,a.qz[i]*a.q[i]]
    fig.add_trace(go.Scatter3d(
        x=x_p,
        y=y_p,
        z=z_p,
        mode='lines',
#        marker_size=1,
        name='pericentro',
        showlegend=False,
#        opacity=0.5,
        marker=dict(
            size=2,
            color='rgb(255, 0, 0)',        # set color to an array/list of desired values
            opacity=0.8
        ),
        line=dict(
            color='rgb(255, 0, 0)',
            width=2
        )
    )
                 )

# Marcador da posição atual
    x_p=[0.0,orb1.x_f[100]]
    y_p=[0.0,orb1.y_f[100]]
    z_p=[0.0,orb1.z_f[100]]
    fig.add_trace(go.Scatter3d(
        x=x_p,
        y=y_p,
        z=z_p,
        mode='lines',
#        marker_size=2,
        name='rnow',
        showlegend=False,
#        opacity=0.5,
        marker=dict(
            size=2,
            color='rgb(0, 255, 0)',        # set color to an array/list of desired values
            opacity=0.8
        ),
        line=dict(
            color='rgb(0, 255, 0)',
            width=2
        )
    )
                 )


# Marcador do vetor velocidade atual
    a=queorb.loc[queorb.index==i]
    x_p=[orb1.x_f[100], orb1.x_f[100] + a.VX_ceu[i]]
    y_p=[orb1.y_f[100], orb1.y_f[100] + a.VY_ceu[i]]
    z_p=[orb1.z_f[100], orb1.z_f[100] + a.VZ_ceu[i]]
    fig.add_trace(go.Scatter3d(
        x=x_p,
        y=y_p,
        z=z_p,
        mode='lines',
#        marker_size=2,
        name='Velocidade',
        showlegend=False,
#        opacity=0.5,
        marker=dict(
            size=2,
            color='rgb(0, 0, 0)',        # set color to an array/list of desired values
            opacity=0.8
        ),
        line=dict(
            color='rgb(0, 0, 0)',
            width=2
        )
    )
                 )


# Marcador do vetor de spin
    a=queorb.loc[queorb.index==i]
#    print(a.spx[i],a.spy[i],a.spz[i])
    x_p=[0.0,a.spx[i]*2*a.q[i]]
    y_p=[0.0,a.spy[i]*2*a.q[i]]
    z_p=[0.0,a.spz[i]*2*a.q[i]]
    fig.add_trace(go.Scatter3d(
        x=x_p,
        y=y_p,
        z=z_p,
        mode='lines',
#        marker_size=1,
        name='spin',
        showlegend=False,
#        opacity=0.5,
        marker=dict(
            size=2,
            color='rgb(255, 0, 0)',        # set color to an array/list of desired values
            opacity=0.8
        ),
        line=dict(
            color='rgb(255, 0, 0)',
            width=2
        )
    )
                 )


# Vetor velocidade no pericentro
    a=queorb.loc[queorb.index==i]
    qq = np.array([a.qx[i], a.qy[i], a.qz[i]])
    ss = np.array([a.spx[i], a.spy[i], a.spz[i]])
    vq = np.cross(ss, qq)
    fig.add_trace(go.Scatter3d(
        x=[a.qx[i]*a.q[i], 10*vq[0] + a.qx[i]*a.q[i]],
        y=[a.qy[i]*a.q[i], 10*vq[1] + a.qy[i]*a.q[i]],
        z=[a.qz[i]*a.q[i], 10*vq[2] + a.qz[i]*a.q[i]],
        mode='lines',
#        marker_size=1,
        name='spin',
        showlegend=False,
#        opacity=0.5,
        marker=dict(
            size=2,
            color='rgb(0, 0, 0)',        # set color to an array/list of desired values
            opacity=0.8
        ),
        line=dict(
            color='rgb(0, 0, 0)',
            width=2
        )
    )
                 )

##################### FIM Loop sobre as órbitas

# Preparo o layout
fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=10, range=[-120,120],),
        yaxis = dict(nticks=10, range=[-120,120],),
        zaxis = dict(nticks=10, range=[-120,120],),
        aspectmode="cube"),
    width=1500,
    margin=dict(r=0, l=50, b=20, t=0),

)


name = 'default'
# Default parameters which are used when `layout.scene.camera` is not provided
camera = dict(
    up=dict(x=0, y=0, z=1),
    center=dict(x=0, y=0, z=0),
    eye=dict(x=1.25, y=-1.25, z=1.25),
    projection=dict(type="orthographic")
)

fig.update_layout(scene_camera=camera, title=name)

fig.show()
