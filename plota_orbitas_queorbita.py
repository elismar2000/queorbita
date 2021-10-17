# Programa para fazer animacão com as órbitas do queorbita
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
from astropy import table
import pandas as pd                 # Para manipular arquivos: read_csv
import plotly.graph_objects as go   # Para plot em 3d
import plotly.express as px
import subprocess                   # Para rodar programas externos
import os                           # Para rodar programas externos
import sys

# Com Pandas, cria o DataFrame queorb
queorb = table.Table.read('/home/elismar/Documentos/Fisica/IC/queorbita/orbits_3rd_attempt/selected_orbits_e1.0.txt', format='ascii')
queorb = queorb.to_pandas()
# Comandos pra ter informações sobre o DataFrame queorb:
#queorb.dtypes
#queorb.head()
#queorb.tail()
#print(queorb)

dirName = "/home/elismar/Documentos/Fisica/IC/queorbita/orbits_3rd_attempt/orbits_temp_e1.0"

# Qual orbita plotar?
#   f=com dyn fric
#   nf=sem dyn fric
#   kepl=kepleriana
#
# plotar=["f","nf","kepl"]   # todas
# plotar=["nf"]
# plotar=["kepl"]
# plotar=["f","nf"]
plotar=["f"]


# Plotar marcadores?
# marcador=["spin","q","Vq","Vnow","Rnow"]   # todas
#marcador=["spin"]                      # só vetor de spin
#marcador=["q"]                         # só vetor pericentro
#marcador=["Vq"]                        # só velocidade no pericentro
#marcador=["Vnow"]                      # só velocidade atual
marcador=["Rnow"]                      # só posição atual



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


if os.path.exists(dirName):
    print("Directory " , dirName ,  " exists.  ")
else:
    sys.exit("Directory" + dirName +  "does not exist. Please run calcula_orbitas_queorbita.py first.")

# Começo a povoar o gráfico com o marcador do centro do sistema de coordenadas:

##############################################################################
# Marcador do centro do sistema de coordenadas (comum a todas as órbitas)
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
##############################################################################


# Escolho as órbitas que serão plotadas, vários exemplos abaixo:
#orb_N = np.arange(1000,1100,1)
#orb_N = [2608,2628,2824,2825,2801]    # No arquivo do queorbita --> numero da linha menos 53
#orb_N = [2799, 2801, 2803, 2805, 2807, 2809, 2811, 2813, 2814, 2816, 2822, 2823, 2824, 2825] # orbitas com e=1 q=2.5

# Abaixo o exemplo de como selecionar órbitas por propriedades, criando um novo dataframe "a"
#a=queorb.loc[(queorb["e"] == 0.9) & (queorb["q"] == 5.7) & (queorb["spy"] >= 0.8) & (queorb["PERIC"] == "POS") & (queorb["Dir"] == "Pro")]
#a=queorb.loc[ (queorb["spy"] >= 0.8) & (queorb["PERIC"] == "POS") & (queorb["Dir"] == "Pro")]


#Para plotar todas as órbitas do DataFrame:
orb_N = queorb.index

# Loop sobre todas as órbitas pra desenhar cada uma delas:
for i in orb_N:
    fname = dirName + "/" + "orb_%05d.dat" % i
    orb = pd.read_csv(fname, sep="\t")
    orb1 = orb                      # Pega somente 1Gyr antes e 1Gyr depois do pericentro (linha 10000).
                                    # Os índices ficam os mesmos do dataframe antigo
    orb1 = orb1.reset_index(drop=True)     # Arruma os índices do orb1 em sequência a partir do 0


# Adiciono colunas no dataframe "orb1" com os modulos do raio R e velocidade V em cada posição :
    orb1["R_f"]=np.sqrt(orb1.x_f**2 + orb1.y_f**2 + orb1.z_f**2)
    orb1["V_f"]=np.sqrt(orb1.Vx_f**2 + orb1.Vy_f**2 + orb1.Vz_f**2)

    orb1["R_nf"]=np.sqrt(orb1.x_nf**2 + orb1.y_nf**2 + orb1.z_nf**2)
    orb1["V_nf"]=np.sqrt(orb1.Vx_nf**2 + orb1.Vy_nf**2 + orb1.Vz_nf**2)

    orb1["R_kepl"]=np.sqrt(orb1.x_kepl**2 + orb1.y_kepl**2 + orb1.z_kepl**2)
    orb1["V_kepl"]=np.sqrt(orb1.Vx_kepl**2 + orb1.Vy_kepl**2 + orb1.Vz_kepl**2)

# Calculo o índice do pericentro em cada tipo de órbita:
    idx_Peri_f=int(orb1[['R_f']].idxmin())
    idx_Peri_nf=int(orb1[['R_nf']].idxmin())
    idx_Peri_kepl=int(orb1[['R_kepl']].idxmin())

# Calculo o índice da posição atual:
    inow=orb1["time"].count() / 2





##############################################################################
# ploto coisas comuns a todas as versões da órbita:
    a=queorb.loc[queorb.index==i]
    # Marcador do vetor de spin
    if "spin" in marcador:
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


    # Marcador da posição atual
    if "Rnow" in marcador:
        x_p=[0.0,a.gx[i]]
        y_p=[0.0,a.gy[i]]
        z_p=[0.0,a.gz[i]]
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
##############################################################################









##############################################################################
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

# Marcador do pericentro
        if "q" in marcador:
        #    print(a.qx[i],a.qy[i],a.qz[i],a.q[i])
            x_p=[0.0, orb1.x_f[idx_Peri_f]]
            y_p=[0.0, orb1.y_f[idx_Peri_f]]
            z_p=[0.0, orb1.z_f[idx_Peri_f]]
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

# Marcador do vetor velocidade atual
        if "Vnow" in marcador:
            x_p=[orb1.x_f[inow], 0.02*orb1.Vx_f[inow] + orb1.x_f[inow]]
            y_p=[orb1.y_f[inow], 0.02*orb1.Vy_f[inow] + orb1.y_f[inow]]
            z_p=[orb1.z_f[inow], 0.02*orb1.Vz_f[inow] + orb1.z_f[inow]]
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
                    color='rgb(200, 0, 0)',        # set color to an array/list of desired values
                    opacity=0.8
                ),
                line=dict(
                    color='rgb(200, 0, 0)',
                    width=4
                )
            )
                         )

# Vetor velocidade no pericentro
        if "Vq" in marcador:
            x_p=[orb1.x_f[idx_Peri_f], 0.02*orb1.Vx_f[idx_Peri_f] + orb1.x_f[idx_Peri_f]]
            y_p=[orb1.y_f[idx_Peri_f], 0.02*orb1.Vy_f[idx_Peri_f] + orb1.y_f[idx_Peri_f]]
            z_p=[orb1.z_f[idx_Peri_f], 0.02*orb1.Vz_f[idx_Peri_f] + orb1.z_f[idx_Peri_f]]
            fig.add_trace(go.Scatter3d(
                x=x_p,
                y=y_p,
                z=z_p,
                mode='lines',
                name='V_q',
                marker=dict(
                    size=2,
                    color='rgb(200, 0, 0)',        # set color to an array/list of desired values
                    opacity=0.8
                ),
                line=dict(
                    color='rgb(200, 0, 0)',
                    width=4
                )
            )
                         )

# FIM Com dynamical friction
##############################################################################









##############################################################################
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



# Marcador do pericentro
        if "q" in marcador:
        #    print(a.qx[i],a.qy[i],a.qz[i],a.q[i])
            x_p=[0.0, orb1.x_nf[idx_Peri_nf]]
            y_p=[0.0, orb1.y_nf[idx_Peri_nf]]
            z_p=[0.0, orb1.z_nf[idx_Peri_nf]]
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

# Marcador do vetor velocidade atual
        if "Vnow" in marcador:
            x_p=[orb1.x_nf[inow], 0.02*orb1.Vx_nf[inow] + orb1.x_nf[inow]]
            y_p=[orb1.y_nf[inow], 0.02*orb1.Vy_nf[inow] + orb1.y_nf[inow]]
            z_p=[orb1.z_nf[inow], 0.02*orb1.Vz_nf[inow] + orb1.z_nf[inow]]
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
                    color='rgb(100, 0, 0)',        # set color to an array/list of desired values
                    opacity=0.8
                ),
                line=dict(
                    color='rgb(100, 0, 0)',
                    width=4
                )
            )
                         )

# Vetor velocidade no pericentro
        if "Vq" in marcador:
            x_p=[orb1.x_nf[idx_Peri_nf], 0.02*orb1.Vx_nf[idx_Peri_nf] + orb1.x_nf[idx_Peri_nf]]
            y_p=[orb1.y_nf[idx_Peri_nf], 0.02*orb1.Vy_nf[idx_Peri_nf] + orb1.y_nf[idx_Peri_nf]]
            z_p=[orb1.z_nf[idx_Peri_nf], 0.02*orb1.Vz_nf[idx_Peri_nf] + orb1.z_nf[idx_Peri_nf]]
            fig.add_trace(go.Scatter3d(
                x=x_p,
                y=y_p,
                z=z_p,
                mode='lines',
                name='V_q',
                marker=dict(
                    size=2,
                    color='rgb(100, 0, 0)',        # set color to an array/list of desired values
                    opacity=0.8
                ),
                line=dict(
                    color='rgb(100, 0, 0)',
                    width=4
                )
            )
                         )

# FIM Sem dynamical friction
##############################################################################










##############################################################################
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
#        color=orb1.time,d
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


# Marcador do pericentro
        if "q" in marcador:
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

# Marcador do vetor velocidade atual
        if "Vnow" in marcador:
            x_p=[orb1.x_kepl[inow], 0.02*orb1.Vx_kepl[inow] + orb1.x_kepl[inow]]
            y_p=[orb1.y_kepl[inow], 0.02*orb1.Vy_kepl[inow] + orb1.y_kepl[inow]]
            z_p=[orb1.z_kepl[inow], 0.02*orb1.Vz_kepl[inow] + orb1.z_kepl[inow]]
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
                    width=4
                )
            )
                         )

# Vetor velocidade no pericentro
        if "Vq" in marcador:
            qq = np.array([a.qx[i], a.qy[i], a.qz[i]])
            ss = np.array([a.spx[i], a.spy[i], a.spz[i]])
            vq = np.cross(ss, qq)
            fig.add_trace(go.Scatter3d(
                x=[a.qx[i]*a.q[i], 0.02*a.vq[i]*vq[0] + a.qx[i]*a.q[i]],
                y=[a.qy[i]*a.q[i], 0.02*a.vq[i]*vq[1] + a.qy[i]*a.q[i]],
                z=[a.qz[i]*a.q[i], 0.02*a.vq[i]*vq[2] + a.qz[i]*a.q[i]],
                mode='lines',
        #        marker_size=1,
                name='V_q',
                showlegend=False,
        #        opacity=0.5,
                marker=dict(
                    size=2,
                    color='rgb(0, 0, 0)',        # set color to an array/list of desired values
                    opacity=0.8
                ),
                line=dict(
                    color='rgb(0, 0, 0)',
                    width=4
                )
            )
                         )




# FIM Kepleriana






















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
