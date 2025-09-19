# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 08:15:52 2025

@author: IsabelCamila
Graphs for subsection All Configurations
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
inicio = 'C:\\Users\\IsabelCamila.LAPTOP-G2J65NSG\\Documents\\disc\\HP\\Python Scripts\\para el git\\'
os.chdir(inicio)

from funcionesECA import Run_all_configurations,Graph_all_configurations, Graph_varianza
from funcionesECA import calcular_mag_ene,obtenerPyDelta
#%% fixed variables
esquemas = ["sequential","block parallel","local clocks", "block sequential"]#,
reglas = [54,90, 150, 110]

#since we are going to need the density and energy of each configuration, we will calculate
#from the start to be efficient
n = 16

#%% run experiment

[mag, ene] = calcular_mag_ene(n)
Run_all_configurations(n, mag, ene, esquemas, B = 3, p = 2, rules = reglas)
#since we want more examples of block sequential and local clocks with different
#number of blocks and maximum periods, we run those again.
Run_all_configurations(n, mag, ene, ["block sequential", "local clocks"], 
                       B = 4, p = 4, rules = reglas)
Run_all_configurations(n, mag, ene, ["block sequential", "local clocks"], 
                       B = 5, p = 5, rules = reglas)


#%% importar datos de density y energy para todas las configuraciones de tamaño n
s = 2**n
df_ene = pd.read_csv("energia.csv", index_col=("Unnamed: 0"), sep = ",")
df_ene = df_ene[df_ene["Tamano Muestra"] == s]
df_ene = df_ene.drop(columns = ["Largo Anillo", "Tamano Muestra"])
df_mag = pd.read_csv("magnetizacion.csv", index_col=("Unnamed: 0"), sep = ",")
df_mag = df_mag[df_mag["Tamano Muestra"] == s]
df_mag = df_mag.drop(columns = ["Largo Anillo", "Tamano Muestra"])

path = inicio

fig_size = [2.5,1.6]
# fig_size = [6,3.6]
#%%
#to graph magnetization and energy for each update scheme, along with each respective rule
for esquema in esquemas:
    Graph_all_configurations(esquema, reglas, df_mag, df_ene, path,
                             limites = [(-.05,.9),(-1.05,.55)], fig_size = fig_size,
                             iteraciones = 150, mostrar = True)
# Graph_varianza(n, path)#, mostrar = True)
#%%
r = 110
scheme = "block sequential"
df_magnetizacion = df_mag
df_energia = df_ene

df_mag_r = df_magnetizacion[df_magnetizacion["regla"] == str(r)+' '+scheme].reset_index().drop(columns = ["regla","index"])
df_ene_r = df_energia[df_energia["regla"] == str(r)+' '+scheme].reset_index().drop(columns = ["regla","index"])
df_mag_0 = {}
df_mag_0[5] = df_mag_r[df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"])
df_mag_0[4] = df_mag_r[~df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_mag_0[4] = df_mag_0[4][df_mag_0[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:]
df_mag_0[3] = df_mag_r[~df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_mag_0[3] = df_mag_0[3][~df_mag_0[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:]

df_ene_0 = {}
df_ene_0[5] = df_ene_r[df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"])
df_ene_0[4] = df_ene_r[~df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_ene_0[4] = df_ene_0[4][df_ene_0[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:]
df_ene_0[3] = df_ene_r[~df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_ene_0[3] = df_ene_0[3][~df_ene_0[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:]

# for B in df_mag_0.keys():
B = 3
fig, ax1 = plt.subplots(figsize=(8,5))#, layout='constrained')
ax2 = ax1.twinx()
ax3 = ax1.twinx()
# ax1.set_ylim((-.1,1.1))
# ax2.set_ylim((-1.1,1.1))
# ax1.set_xlabel("Frecuency")
# ax1.set_ylabel("Density")
# ax2.set_ylabel("Energy")
xx = [int(i) for i in df_mag_0[B].columns]
p1 = []
p2 = []
p3 = []
for i in range(len(df_mag_0[B])):
# i = 0
    p1.append(np.mean(df_mag_0[3].iloc[i]))
    p2.append(np.mean(df_mag_0[4].iloc[i]))
    p3.append(np.mean(df_mag_0[5].iloc[i]))
    
counts, bins = np.histogram(p1)
ax1.stairs(counts,bins, orientation = "vertical", fill = True, label = "3 blocks")
counts, bins = np.histogram(p2)
ax1.stairs(counts, bins, orientation = "vertical", ls = ":", color = 'k', fill = True, alpha = .5, label = "4 blocks")
counts, bins = np.histogram(p3)
ax1.stairs(counts, bins, orientation = "vertical", ls = ":", color = 'm', fill = True, alpha = .5, label = "5 blocks")
plt.title("Density\nRule "+str(r)+' '+scheme+' ring '+str(n)+"\n Average over different schemes with "+str(B)+" blocks")
ax1.legend(loc="lower right", bbox_to_anchor=(.75, -.05), 
           bbox_transform=fig.transFigure, ncol=3)
# ax2.legend(loc="lower right", bbox_transform=fig.transFigure, ncol=3)
# plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'mag'+str(B)+'.pdf', bbox_inches='tight')
# plt.draw()
# if mostrar:
plt.show()
plt.close()


#%%Limpiar Datos
df_ene = pd.read_csv("energia - copia.csv", index_col=("Unnamed: 0"), sep = ",")
df_mag = pd.read_csv("magnetizacion - copia.csv", index_col=("Unnamed: 0"), sep = ",")
df_mag = df_mag[df_mag["Tamano Muestra"] == 2**16].reset_index().drop(columns = ["index"])
df_ene = df_ene[df_ene["Tamano Muestra"] == 2**16].reset_index().drop(columns = ["index"])
#%%block sequential
df_mag_r = df_mag[df_mag["regla"].str.contains("block sequential")].reset_index().drop(columns = ["index"])
df_ene_r = df_ene[df_ene["regla"].str.contains("block sequential")].reset_index().drop(columns = ["index"])

df_mag_0 = {}
df_mag_0[5] = df_mag_r[df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_mag_0[4] = df_mag_r[~df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_mag_0[4] = df_mag_0[4][df_mag_0[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])[-128:]
df_mag_0[3] = df_mag_r[~df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_mag_0[3] = df_mag_0[3][~df_mag_0[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])[-128:]
df_ene_0 = {}
df_ene_0[5] = df_ene_r[df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_ene_0[4] = df_ene_r[~df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_ene_0[4] = df_ene_0[4][df_ene_0[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])[-128:]
df_ene_0[3] = df_ene_r[~df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
df_ene_0[3] = df_ene_0[3][~df_ene_0[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])[-128:]

#%% block parallel

df_mag_BP = df_mag[df_mag["regla"].str.contains("block parallel")].reset_index().drop(columns = ["index"])
df_mag_BP = df_mag_BP.reset_index().drop(columns = ["index"]).iloc[-128:].reset_index().drop(columns = ["index"])
df_ene_BP = df_ene[df_ene["regla"].str.contains("block parallel")].reset_index().drop(columns = ["index"])
df_ene_BP = df_ene_BP.reset_index().drop(columns = ["index"]).iloc[-128:].reset_index().drop(columns = ["index"])
#%% sequential
df_mag_SEQ = df_mag[df_mag["regla"].str.contains("sequential")].reset_index().drop(columns = ["index"])
df_mag_SEQ = df_mag_SEQ[~df_mag_SEQ["regla"].str.contains("c")].reset_index().drop(columns = ["index"])
df_mag_SEQ = df_mag_SEQ.reset_index().drop(columns = ["index"]).iloc[-128:].reset_index().drop(columns = ["index"])
df_ene_SEQ = df_ene[df_ene["regla"].str.contains("sequential")].reset_index().drop(columns = ["index"])
df_ene_SEQ = df_ene_SEQ[~df_ene_SEQ["regla"].str.contains("c")].reset_index().drop(columns = ["index"])
df_ene_SEQ = df_ene_SEQ.reset_index().drop(columns = ["index"]).iloc[-128:].reset_index().drop(columns = ["index"])
#%% local clocks
df_mag_r = obtenerPyDelta(df_mag,16)
df_ene_r = obtenerPyDelta(df_ene,16)
df_mag_LC = df_mag_r[df_mag_r["regla"].str.contains("local clocks")].reset_index().drop(columns = ["index"])
df_ene_LC = df_ene_r[df_ene_r["regla"].str.contains("local clocks")].reset_index().drop(columns = ["index"])
esquemas = [elemento for elemento in list(df_mag_LC["Esquema"])]
periodos = {str(esquema): max(esquema[1][0]) for esquema in esquemas}
dic_mag = {}
dic_ene = {}
for periodo in [2,4,5]:
    aux = [esquema for esquema in esquemas if periodos[str(esquema)] == periodo]
    dic_mag[periodo] = df_mag_LC[df_mag_LC["Esquema"].isin(aux)].reset_index().drop(columns = ["index"])[:128]
    dic_ene[periodo] = df_ene_LC[df_ene_LC["Esquema"].isin(aux)].reset_index().drop(columns = ["index"])[:128]

#%% juntar todo

df_ene = pd.concat([df_ene_0[3],df_ene_0[4]], ignore_index= True)
df_ene = pd.concat([df_ene,df_ene_0[5]], ignore_index= True)
df_ene = pd.concat([df_ene,df_ene_BP], ignore_index= True)
df_ene = pd.concat([df_ene,df_ene_SEQ], ignore_index= True)
df_ene = pd.concat([df_ene,dic_ene[2]], ignore_index= True)
df_ene = pd.concat([df_ene,dic_ene[4]], ignore_index= True)
df_ene = pd.concat([df_ene,dic_ene[5]], ignore_index= True)

df_mag = pd.concat([df_mag_0[3],df_mag_0[4]], ignore_index= True)
df_mag = pd.concat([df_mag,df_mag_0[5]], ignore_index= True)
df_mag = pd.concat([df_mag,df_mag_BP], ignore_index= True)
df_mag = pd.concat([df_mag,df_mag_SEQ], ignore_index= True)
df_mag = pd.concat([df_mag,dic_mag[2]], ignore_index= True)
df_mag = pd.concat([df_mag,dic_mag[4]], ignore_index= True)
df_mag = pd.concat([df_mag,dic_mag[5]], ignore_index= True)
#%%guardar
df_ene.to_csv("energia.csv")
df_mag.to_csv("magnetizacion.csv")


