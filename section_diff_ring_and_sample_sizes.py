# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:06:02 2025

@author: IsabelCamila
Graphs for subsection Different ring and sample sizes
"""

import pandas as pd
import os
#please write the location of the files here
inicio = 'C:\\Users\\IsabelCamila.LAPTOP-G2J65NSG\\Documents\\disc\\HP\\Python Scripts\\para el git\\'
os.chdir(inicio)
from funcionesECA import bin_r, run_experiments_samples, Graph_Ring_and_Sample_Sizes
reglas = {regla:bin_r(regla) for regla in [54,110,150,90]}
esquemas = ["sequential", "block parallel","block sequential","local clocks"]#,
#%% In this section we execute the experiments
#first we must retrive the sample
n_muestra = 0
df_muestra = pd.read_csv("muestra_32_8.csv", index_col=("Unnamed: 0"))
muestra_32_8 = [int(i) for i in df_muestra[str(n_muestra)]]
df_muestra = pd.read_csv("muestra_128_8.csv", index_col=("Unnamed: 0"))
muestra_128_8 = [int(i) for i in df_muestra[str(n_muestra)]]
df_muestra = pd.read_csv("muestra_32_38.csv", index_col=("Unnamed: 0"))
muestra_32_38 = [int(i) for i in df_muestra[str(n_muestra)]]
df_muestra = pd.read_csv("muestra_128_38.csv", index_col=("Unnamed: 0"))
muestra_128_38 = [int(i) for i in df_muestra[str(n_muestra)]]
df_muestra = pd.read_csv("muestra_32_138.csv", index_col=("Unnamed: 0"))
muestra_32_138 = [int(i) for i in df_muestra[str(n_muestra)]]
df_muestra = pd.read_csv("muestra_128_138.csv", index_col=("Unnamed: 0"))
muestra_128_138 = [int(i) for i in df_muestra[str(n_muestra)]]

#%% experiments with a given sample
run_experiments_samples(n = 8, reglas = reglas, muestra = muestra_32_8, esquemas = esquemas)
run_experiments_samples(n = 8, reglas = reglas, muestra = muestra_128_8, esquemas = esquemas)
run_experiments_samples(n = 38, reglas = reglas, muestra = muestra_32_38, esquemas = esquemas)
run_experiments_samples(n = 38, reglas = reglas, muestra = muestra_128_38, esquemas = esquemas)
run_experiments_samples(n = 138, reglas = reglas, muestra = muestra_32_138, esquemas = esquemas)
run_experiments_samples(n = 138, reglas = reglas, muestra = muestra_128_138, esquemas = esquemas)

#%% In this section we create the graphs

df_mag = {}
df_ene = {}
for n in [8,38,138]:
    df_mag[n] = pd.read_csv("magnetizacion_muestras_"+str(n)+".csv", index_col=("Unnamed: 0"), sep = ",")
    df_ene[n] = pd.read_csv("energia_muestras_"+str(n)+".csv", index_col=("Unnamed: 0"), sep = ",")

path = inicio
fig_size = [5.5,3.3]

Graph_Ring_and_Sample_Sizes([df_mag,df_ene], rule = 54, list_ylimits = [(0,.9),(-1,.8)], 
                            fig_size = fig_size,
                            camino = path, iteraciones = 150, schemes = esquemas, to_show = True)
Graph_Ring_and_Sample_Sizes([df_mag,df_ene], rule = 110, list_ylimits = [(0,.9),(-1,.8)], 
                            fig_size = fig_size,
                            camino = path, iteraciones = 150, schemes = esquemas, to_show = True)

Graph_Ring_and_Sample_Sizes([df_mag,df_ene], rule = 90, list_ylimits =  [(0,.9),(-1,.8)], 
                            fig_size = fig_size,
                            camino = path, iteraciones = 150, schemes = esquemas)#, to_show = True)

Graph_Ring_and_Sample_Sizes([df_mag,df_ene], rule = 150, list_ylimits =  [(0,.9),(-1,.8)], 
                            fig_size = fig_size,
                            camino = path, iteraciones = 150, schemes = esquemas)#, to_show = True)
