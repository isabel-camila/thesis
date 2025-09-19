# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 12:09:01 2025

@author: IsabelCamila

Definition of Functions
"""
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import pandas as pd
from math import lcm,ceil

def esquema_bs(n,B):
    #n is the size of the ring
    #B is the number of blocks
    dic_aux = {i:[] for i in range(B)} #i create an auxiliar dictionary
    opciones = list(range(n)) #then i list the options
    while len(opciones) != 0:
        b = rd.choice(opciones) #randomly select a cell 
        dic_aux[rd.choice(range(B))].append(b) #assign it to a block
        opciones.remove(b) #prevent repetition
    dic_aux = {i:dic_aux[i] for i in dic_aux.keys() if len(dic_aux[i])!=0} #remove possible empty blocks
    return [{i:dic_aux[list(dic_aux.keys())[i]] for i in range(len(dic_aux))},
            {i:dic_aux[list(dic_aux.keys())[i]] for i in range(len(dic_aux))}] 
    #return both the "definition" of the update mode and the associated blocks

def esquema_bp(n,B):
    definicion = esquema_bs(n,B) #we generate B blocks using block seq strategy
    lista = [len(definicion[0][b]) for b in definicion[0].keys()]
    mcm = lcm(*lista) #find the least common multiple
    while mcm > 1000: #do it as many times as necessary to prevent a number of sub steps too big
        definicion = esquema_bs(n,B) 
        lista = [len(definicion[0][b]) for b in definicion[0].keys()]
        mcm = lcm(*lista)
    bloques = {i:[] for i in range(mcm)}
        
    #now we obtain the sequence of blocks associated to the definition
    for b in range(mcm):
        lista = []
        for celdas in definicion[0].values():
            lista.append([celda for celda in celdas if celdas.index(celda) == b%len(celdas)])
        bloques[b] = [i for sl in lista for i in sl]
    return [bloques,definicion[0]] 
    #note that the sequence of blocks is different from the definition
    #because of Block-Parallel is the dual to Block-Sequential.

def esquema_bp16(n = 16):
    #specifically for ring of size 16, we defined particular sizes of the blocks 
    #of the partition, to ensure that the cycles weren't too too long.
    dic_aux = {i:[] for i in range(4)}
    opciones = list(range(n))
    b = rd.choice(opciones)
    dic_aux[0].append(b)
    opciones.remove(b)
    for i in range(3):
        b = rd.choice(opciones)
        dic_aux[1].append(b)
        opciones.remove(b)
    for i in range(5):
        b = rd.choice(opciones)
        dic_aux[2].append(b)
        opciones.remove(b)
    dic_aux[3] = opciones #generamos los B bloques con bs
    lista = [len(dic_aux[b]) for b in dic_aux.keys()]
    mcm = lcm(*lista)
    aux = {i:[] for i in range(mcm)}
    #procedemos a reorganizar para que se ejecuten en modo bloque paralelo
    for b in range(mcm):
        lista = []
        for celdas in dic_aux.values():
            lista.append([celda for celda in celdas if celdas.index(celda) == b%len(celdas)])
        aux[b] = [i for sl in lista for i in sl]
    return [aux,dic_aux]

def esquema_lc(n,p):
    #n is the size of the ring
    #p is the max period that a cell can have
    relojes = [rd.choice(range(1,p+1)) for i in range(n)] #we create a list of periods
    fases = [rd.choice(range(r)) for r in relojes] #and a list of shifts
    mcm = lcm(*relojes)
    bloques = {i:[] for i in range(mcm)}
    for b in range(mcm): #and now i need the blocks associated with the definition
        bloques[b] = [i for i in range(n) if b%relojes[i]==(fases[i])%relojes[i]]
    return [bloques,[relojes,fases]]

def crear_bloques(n_esquema,n_bloques,periodo,size_ring, especial = False):
    #one function to create the update mode by family
    if n_esquema == "parallel": #all at once
        return  [{0:list(range(size_ring))},{0:list(range(size_ring))}]
    elif n_esquema == "block parallel":
        if especial: #to protect againt too many substeps. only works for n = 16
            return esquema_bp16(size_ring) 
        else: #for the rest.
            return esquema_bp(size_ring, n_bloques)
    elif n_esquema == "sequential":
        lista_aux = [i for i in range(size_ring)]
        rd.shuffle(lista_aux)
        return [{i: [lista_aux[i]] for i in range(size_ring)},{i: [lista_aux[i]] for i in range(size_ring)}]
    elif n_esquema == "sequential c":
        #we called it canonical sequential. it means that the cells update one-by-one 
        #from left to right
        return [{i: [i] for i in range(size_ring)},{i: [i] for i in range(size_ring)}]   
    elif n_esquema == "block sequential":
        return esquema_bs(size_ring, n_bloques)
    elif n_esquema == "local clocks":
        return esquema_lc(size_ring, periodo)
        #bloques = {0: [1, 2, 3, 4, 5, 6, 7, 11, 13, 16, 18, 19], 1: [0, 3, 5, 8, 10, 12, 13, 14, 16, 18, 19], 2: [1, 2, 3, 4, 5, 9, 11, 13, 15, 16, 18, 19], 3: [3, 5, 6, 10, 13, 14, 16, 17, 18, 19], 4: [0, 1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 16, 18, 19], 5: [3, 5, 9, 10, 13, 14, 16, 18, 19], 6: [1, 2, 3, 4, 5, 6, 11, 13, 15, 16, 18, 19], 7: [0, 3, 5, 8, 10, 12, 13, 14, 16, 17, 18, 19], 8: [1, 2, 3, 4, 5, 7, 9, 11, 13, 16, 18, 19], 9: [3, 5, 6, 10, 13, 14, 16, 18, 19], 10: [0, 1, 2, 3, 4, 5, 8, 11, 12, 13, 15, 16, 18, 19], 11: [3, 5, 9, 10, 13, 14, 16, 17, 18, 19]}
    elif n_esquema == "bipartite":
        return [{0: range(0,size_ring,2), 1: range(1,size_ring,2)},{0: range(0,size_ring,2), 1: range(1,size_ring,2)}]

def bin_r(rule):
    #takes a rule in base 10 and translates it to base2, with which we obtain the
    #definition of the rule
    rule_b = '0'*8+ bin(rule)[2:]
    rule_b = rule_b[len(rule_b)-8:]
    return rule_b

def bin_c(conf,n):
    #takes the "name" of a configuration in base10 and translates it to base2
    #that way i can waste less storage memory
    aux0 = bin(conf)[2:]
    if len(aux0) < n:
        aux0 = '0'*(n-len(aux0))+aux0
    return aux0
def apply_rule(conf, rule_b, n, Bloque):
    #we apply the rule one time only on the cells that belong to the Bloque
    #conf must be received in base10
    aux0 = bin_c(conf, n) #then we change it to base2
    aux1 = ''
    for i in range(n):
        if i in Bloque: #here we ensure that we respect the block
            V = int(aux0[i-1]+aux0[i]+aux0[(i+1)%n],2)
            aux1 = aux1 + rule_b[len(rule_b)-V-1]
        else:
            aux1 = aux1 + aux0[i]
    return int(aux1,2)

def ejecutar_dinamicas(rule,nn,bloquess,itera,conf):
    #rule is the ECA rule in binary
    #nn is the length of the ring. must be a natural number
    #bloquess is a dictionary in which
        #the keys are the substeps
        #the values are the cells in the block of the substep
    #itera number of iterations
    #conf is the configuration to be updated, in base 10
    dinamicas = {}
    all_substeps = [conf]
    full_steps = [conf]
    for iteracion in range(itera):
        for subpaso in bloquess.keys():
            #we find the image after each substep
            all_substeps.append(apply_rule(all_substeps[-1], rule, nn, bloquess[subpaso]))   
        #we find the final step after all the substeps are done
        full_steps.append(all_substeps[-1])
        if (all_substeps[-1] in full_steps[:-1]): #we check if the newest image has already been found (ie, we have a cycle)
            #this is a pointer, to find where the cycle begins
            point = [i for i,x in enumerate(full_steps) if x == all_substeps[-1]]
            #this saves the cycle
            dinamicas["ciclo"] = full_steps[point[0]:]
            #we stop with the iterations and instead we repeat the cycle as many times as necessary
            full_steps = full_steps + ceil(itera/len(dinamicas["ciclo"]))*dinamicas["ciclo"] 
            break    #once we find a cycle, we get out of the loop.
    #in case the total number of steps ends up being longer thanthe established number of iterations
    #this can happen if itera/len(dinamicas["ciclo"])*dinamicas["ciclo"] is longer than itera
    dinamicas["full steps"] = full_steps[:itera+1]
    dinamicas["bloques"] = bloquess
    return dinamicas

def ver_dinamica(din,n, ver = True):
    #din must be a list with the dynamic
    #ver can be False or True. it determines whether the din is visualized
    #this function creates an array that contains y the time component and x the configuration.
    visualizacion = []
    for elemento in din:
        visualizacion.append([int(i) for i in list(bin_c(elemento,n))])
    visualizacion = np.array(visualizacion)
    
    if ver:
        plt.matshow(visualizacion, cmap = 'Greens')
    #cmap can be: 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
    #                  'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
    #                  'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
    #                  'viridis', 'plasma', 'inferno', 'magma', 'cividis'
    return visualizacion
def calcular_magnetizacion(visualizacion, nn, ver = True):
    #we calculate the magnetization of a given dinamics, through the array obtained after using ver_dinamica.
    magnets = np.sum(visualizacion, axis=1)/nn
    if ver:
        plt.plot(range(len(magnets)),magnets)
        plt.title("Magnetization",)
        plt.xlabel("Number of Iteration")
        plt.ylabel("Amount of Ones")
        plt.show()
    return magnets
def avg_magnetizacion(dict_dinamica, titulo, nn, t_muestra, itera, ver = True):
    densidades = np.array([])
    for configuracion in dict_dinamica:
        vis = ver_dinamica(dict_dinamica[configuracion]["full steps"], nn, ver = False)
        densidad = calcular_magnetizacion(vis, nn, ver = False)
        densidades = np.concatenate([densidades,densidad])
    densidades = np.reshape(densidades, (len(dict_dinamica),itera+1)).T
    avg_densidad = np.mean(densidades, axis = 1)
    avg_densidad = np.round(avg_densidad,2)
    if ver:
        plt.plot(range(len(avg_densidad)),avg_densidad)
        plt.title("Magnetization\n" + titulo+ " ring "+ str(nn)+ "\n sample size "+str(t_muestra))
        plt.xlabel("Number of Iteration")
        plt.ylabel("Density")
        plt.show()
    return avg_densidad

def calcular_energia(visualizacion, nn, ver = True):
    #recibe un array con las dinámicas cada iteración.
    #dibuja una progresión de la energía según se aplica la regla.
    lista = []
    nn = len(visualizacion[0])
    for t in range(len(visualizacion)):
        Suma = 0
        for i in range(len(visualizacion[t])):
            suma = -(-1+2*visualizacion[t][i])*((-1+2*visualizacion[t][i-1])+(-1+2*visualizacion[t][(i+1)%nn]))/2
            Suma += suma
        lista.append(Suma)
    if ver:
        plt.plot(range(len(visualizacion)), lista, 'o')
        plt.plot(range(len(visualizacion)), lista)
        plt.title("Energy",)
        plt.xlabel("Number of Iteration")
        plt.show()
        # plt.savefig('ene.png')
    return lista
def avg_energia(dict_dinamica, titulo, nn, t_muestra, itera, ver = True):
    energias = np.array([])
    for configuracion in dict_dinamica:
        vis = ver_dinamica(dict_dinamica[configuracion]["full steps"], nn, ver = False)
        energy = calcular_energia(vis, nn, ver = False)
        energias = np.concatenate((energias,energy))
    energias = np.reshape(energias, (len(dict_dinamica),itera+1)).T
    avg_energy = np.mean(energias, axis = 1)
    avg_energy = np.round(avg_energy,2)
    if ver:
        plt.plot(range(len(avg_energy)),avg_energy)
        plt.title("Energy\n" + titulo + " ring "+ str(nn) + "\n sample size "+str(t_muestra))
        plt.xlabel("Number of Iteration")
        plt.show()
    return avg_energy
def run_experiments_samples(n, reglas, muestra, esquemas, n_bloques = 3, periodo = 5):
    try: 
        df_mag = pd.read_csv("magnetizacion_muestras_"+str(n)+".csv", index_col=("Unnamed: 0"), sep = ",")
    except FileNotFoundError: 
        #If the file of magnetization for a ring of size n is not found it creates a new one
        df_mag = pd.DataFrame()
    try: 
        df_ene = pd.read_csv("energia_muestras_"+str(n)+".csv", index_col=("Unnamed: 0"), sep = ",")
    except FileNotFoundError: #If the file ofenergy does not exist in the current location
        df_ene = pd.DataFrame()
    iteraciones = 1000
    #we perform 32 iterations for each update mode, for each rule
    for i in range(32):
        #we will perform the experiment with 32 different update modes from each family
        dinamica = {}
        promedio_mag = {}
        promedio_ene = {}
        ejecutados = []
        bloques = {}
        for esquema in esquemas:
            # we need to define the function for each update mode
            bloques[esquema] = crear_bloques(esquema, n_bloques = 3, periodo = 5, size_ring = n)    
        for esquema in esquemas:
            for regla in reglas.keys():
                llave = str(regla)+' '+esquema
                dinamica[llave] = {}
                for conf in muestra:
                    #aquí separé el proceso para cuando tenga que keyboard interrupt no pierda todo mi progreso
                    dinamica[llave][conf] = ejecutar_dinamicas(reglas[regla], n, bloques[esquema][0], iteraciones, conf) 
                promedio_mag[llave] = avg_magnetizacion(dinamica[llave], llave,n, len(muestra), iteraciones, ver = False)
                promedio_ene[llave] = avg_energia(dinamica[llave], llave, n, len(muestra), iteraciones, ver = False)
                
            ejecutados.append(esquema)
        
        df = pd.DataFrame(promedio_mag, columns = promedio_mag.keys()).T
        df = df.reset_index()
        df = df.rename(columns = {"index":"regla"})
        df = df.drop(range(len(reglas)*len(ejecutados),len(df)))
        df["Largo Anillo"] = [n for i in range(len(reglas)*len(ejecutados))]
        df["Tamano Muestra"] = [len(muestra) for i in range(len(reglas)*len(ejecutados))]
        df["Esquema"] = [[esquema,bloques[esquema][1]] for esquema in ejecutados for regla in reglas]
        #we ensure that the columns are in the right order
        cols = df.columns.tolist()
        cols = cols[-4:]+cols[:-4]
        df = df[cols]
        
        df_mag = pd.concat([df_mag,df],ignore_index=True)
        
        df = pd.DataFrame(promedio_ene, columns = promedio_ene.keys()).T
        df = df.reset_index()
        df = df.rename(columns = {"index":"regla"})
        df = df.drop(range(len(reglas)*len(ejecutados),len(df)))
        df["Largo Anillo"] = [n for i in range(len(reglas)*len(ejecutados))]
        df["Tamano Muestra"] = [len(muestra) for i in range(len(reglas)*len(ejecutados))]
        df["Esquema"] = [[esquema,bloques[esquema][1]] for esquema in ejecutados for regla in reglas]
        cols = df.columns.tolist()
        cols = cols[-4:]+cols[:-4]
        df = df[cols]
        
        df_ene = pd.concat([df_ene,df],ignore_index=True)
    df_mag.to_csv("magnetizacion_muestras_"+str(n)+".csv") #save the data
    df_ene.to_csv("energia_muestras_"+str(n)+".csv") #save the data

def calcular_mag_ene(nn):
    #calculamos mag y ene para todas las configuraciones de un anillo de tamaño nn.
    #retorna una dos listas, una con density y la otra con energy.
    magnet = {}
    energy = {}
    for conf in range(2**nn):
        conf_b = [int(i) for i in bin_c(conf, nn)]
        mag_conf = np.sum(conf_b)/nn
        ene_conf = 0
        for i in range(nn):
            ene_conf = ene_conf + -(-1+2*conf_b[i])*((-1+2*conf_b[i-1])+(-1+2*conf_b[(i+1)%nn]))/2
        magnet[conf] = np.round(mag_conf,2)
        energy[conf] = np.round(ene_conf,2)
    return [magnet, energy]

def ejecutar_todo_m_y_e(rule, nn, bloquess, mag_dic, ene_dic, itera = 1000):
    #rule es la regla ACE en formato binario
    #nn es el largo del anillo. es un número natural
    #bloquess es un diccionario con el orden en que 
        #las keys son los subpaso
        #los values son las celdas que se actualizan en ese subpaso
    #itera es el número de iteraciones máximo (natural)
    #conf es la configuración que vamos a actualizar, en formato decimal.
    #mag y ene tienen que ser calculador previamente para optimizar tiempo.
    all_substeps = {}
    full_steps = {}
    for conf in range(2**nn):
        all_substeps[conf] = [conf]
        for subpaso in bloquess.keys():
            all_substeps[conf].append(apply_rule(all_substeps[conf][-1], rule, nn, bloquess[subpaso]))
        full_steps[conf] = all_substeps[conf][-1]
    dinamicas = {}
    imanes = np.zeros((itera+1,))
    energias = np.zeros((itera+1,))
    for conf in range(2**nn):
        dinamicas[conf] = {}
        dyns = [conf]
        for iteracion in range(itera):
            dyns.append(full_steps[dyns[-1]])
            if (full_steps[dyns[-1]] in dyns):
                point = [i for i,x in enumerate(dyns) if x == full_steps[dyns[-1]]]
                dinamicas[conf]["transiente"] = dyns[:point[0]]
                try: 
                    if len(dyns[point[0]:])==2:
                        if dyns[point[0]:][0]==dyns[point[0]:][1]:
                            dinamicas[conf]["ciclo"] = dyns[0]
                        else:
                            dinamicas[conf]["ciclo"] = dyns[point[0]:]
                    else:
                        dinamicas[conf]["ciclo"] = dyns[point[0]:]
                except:
                    print(dyns[point[0]:])
                #dejo de calcular la siguiente iteración instead repito el ciclo el número necesario de veces
                #dyns = dyns + ceil(itera/len(dinamicas[conf]["ciclo"]))*dinamicas[conf]["ciclo"] 
                break    #break permite que las iteraciones se detengan al encontrar un ciclo
        
        # dinamicas[conf]["full steps"] = dyns[:itera+1]
        #aquí en lugar de almacenar todas las dinámica de ene y den, voy sumando
        #de esa manera se gasta menos memoria
        imanes = imanes +  np.array([mag_dic[i] for i in dyns[:itera+1]])
        energias = energias + np.array([ene_dic[i] for i in dyns[:itera+1]])
    #para exportar los valores de magnetizacion y energía se van dentro del dict de dinamicas.
    dinamicas["magnet"] = np.round(imanes/(2**nn),2)
    dinamicas["energy"] =  np.round(energias/(2**nn),2)
    
    return dinamicas
def Run_all_configurations(n, mag, ene, esquemas, B, p, rules):
    iteraciones = 1000
    try :
        df_mag = pd.read_csv("magnetizacion.csv", index_col=("Unnamed: 0"))
        df_mag = df_mag.rename(columns = {str(i):int(i) for i in range(iteraciones+1)})
    except FileNotFoundError: 
        #If the file of magnetization does not exist yet
        df_mag = pd.DataFrame()
    try :
        df_ene = pd.read_csv("energia.csv", index_col=("Unnamed: 0"))
        df_ene = df_ene.rename(columns = {str(i):int(i) for i in range(iteraciones+1)})
    except FileNotFoundError: 
        #If the file of energy does not exist yet
        df_mag = pd.DataFrame()
    for i in range(32):
        dinamica = {}
        promedio_mag = {}
        promedio_ene = {}
        ejecutados = []
        bloques = {}
        for esquema in esquemas:
            #the property "especial" forces the length of the blocks to be manageble for
            #the block parallel update mode.
            bloques[esquema] = crear_bloques(esquema, B, p, n, especial = True)
            for rule in rules.keys():
                llave = str(rule)+' '+esquema
                dinamica[llave] = ejecutar_todo_m_y_e(rules[rule], n, bloques[esquema][0], mag, ene, iteraciones) 
                promedio_mag[llave] = dinamica[llave]["magnet"]
                promedio_ene[llave] = dinamica[llave]["energy"]
                
            ejecutados.append(esquema)
        
        df = pd.DataFrame(promedio_mag, columns = promedio_mag.keys()).T
        df = df.reset_index()
        df = df.rename(columns = {"index":"regla"})
        df = df.drop(range(len(rules)*len(ejecutados),len(df)))
        df["Tamano Muestra"] = [2**n for i in range(len(rules)*len(ejecutados))]
        df["Largo Anillo"] = [n for i in range(len(rules*len(ejecutados)))]
        df["Esquema"] = [[esquema,bloques[esquema][1]] for esquema in ejecutados for regla in rules]
        cols = df.columns.tolist()
        cols = cols[-3:]+cols[:-3] #proper order of columns
        df = df[cols]
        df_mag = pd.concat([df_mag,df],ignore_index=True)
        df_mag.to_csv("magnetizacion.csv") #save the data
        
        df = pd.DataFrame(promedio_ene, columns = promedio_ene.keys()).T
        df = df.reset_index()
        df = df.rename(columns = {"index":"regla"})
        df = df.drop(range(len(rules)*len(ejecutados),len(df)))
        df["Tamano Muestra"] = [2**n for i in range(len(rules*len(ejecutados)))]
        df["Largo Anillo"] = [n for i in range(len(rules*len(ejecutados)))]
        df["Esquema"] = [[esquema,bloques[esquema][1]] for esquema in ejecutados for regla in rules]
        cols = df.columns.tolist()
        cols = cols[-3:]+cols[:-3]
        df = df[cols]
        df_ene = pd.concat([df_ene,df],ignore_index=True)
        df_ene.to_csv("energia.csv") #save the data

def Graph_Ring_and_Sample_Sizes(list_df, rule, list_ylimits, camino, fig_size, schemes, iteraciones = 1000, to_show = False):
    path = camino + 'sample\\'
    # list_ylimits = [(0,1),(-1,1)]
    df_mag = list_df[0].copy()
    df_ene = list_df[1].copy()
    colores = {str([8,32]): 'midnightblue', str([38,32]): 'green', str([138,32]): 'maroon',
               str([8,128]): 'lightblue', str([38,128]): 'lightgreen', str([138,128]): 'red'}
    for scheme in schemes:
        df_mag_ = pd.DataFrame()
        df_ene_ = pd.DataFrame()
        for n in [8,38,138]:
            df_aux = df_mag[n][df_mag[n]["regla"] == str(rule)+' '+scheme]
            columnas = ["Largo Anillo", "Tamano Muestra", "Numero de Muestra","Esquema", "regla"] + [str(i) for i in range(iteraciones + 1)]
            df_aux = df_aux[columnas]
            df_aux = df_aux.reset_index().drop(columns = ["index","Numero de Muestra","Esquema","regla"])
            df_mag_ = pd.concat([df_mag_,df_aux], ignore_index=True)
            df_aux = df_ene[n][df_ene[n]["regla"] == str(rule)+' '+scheme]
            columnas = ["Largo Anillo", "Tamano Muestra", "Numero de Muestra","Esquema", "regla"] + [str(i) for i in range(iteraciones + 1)]
            df_aux = df_aux[columnas]
            df_aux = df_aux.reset_index().drop(columns = ["index","Numero de Muestra","Esquema","regla"])
            df_ene_ = pd.concat([df_ene_,df_aux], ignore_index=True)
    
        fig, ax1 = plt.subplots()#figsize=fig_size)#, layout='constrained')
        ax2 = ax1.twinx()
        ax1.set_ylim(list_ylimits[0])
        ax1.set_ylabel("Density", fontsize=9)
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_ylabel("Energy", fontsize=9)
        ax2.set_ylim(list_ylimits[1])
        ax2.tick_params(axis='both', which='major', labelsize=8)
        for n in [8,38,138]:
            for tamano_muestra in [32,128]:
                df_aux = df_mag_[df_mag_["Tamano Muestra"] == tamano_muestra]
                aux = np.mean(df_aux[df_aux["Largo Anillo"] == n].reset_index().drop(columns = ["Tamano Muestra", "Largo Anillo","index"]).iloc[:32].to_numpy(),axis = 0)
                ax1.plot(range(len(aux)), aux, label = 'Ring:'+str(n)+" Sample size:"+str(tamano_muestra), color = colores[str([n,tamano_muestra])], linestyle='-')
                df_aux = df_ene_[df_ene_["Tamano Muestra"] == tamano_muestra]
                aux = np.mean(df_aux[df_aux["Largo Anillo"] == n].reset_index().drop(columns = ["Tamano Muestra", "Largo Anillo","index"]).iloc[:32].to_numpy(),axis = 0)/n
                ax2.plot(range(len(aux)), aux, label = 'Ring:'+str(n)+" Sample size:"+str(tamano_muestra), color = colores[str([n,tamano_muestra])], linestyle='dotted')
                
        plt.title(#"Density and Energy\n "+
                  "Rule "+str(rule)+' under '+scheme.title()+" Update Mode")#"\n Average over different sample size.")
        ax1.legend(title = "Density",loc="center", bbox_to_anchor=(.52, -0.04), 
                    fontsize = 8,
                    bbox_transform=fig.transFigure, columnspacing = 1,handlelength = 1.5,
                    ncol=3).set_visible(False)
        ax2.legend(title = "Energy",loc="center", bbox_to_anchor=(.52, -.23), 
                    fontsize = 8,
                    bbox_transform=fig.transFigure, columnspacing = 1,handlelength = 1.5,
                    ncol=3).set_visible(False)
        ax1.grid(axis = 'y',color = 'gray', linestyle = '--', linewidth = 0.1)
        ax2.grid(axis = 'y',color = 'gray', linestyle = 'dotted', linewidth = 0.1)
        ax1.grid(axis = 'x',color = 'gray', linestyle = '--', linewidth = 0.1)
        fig.set_size_inches(fig_size[0],fig_size[1])
        plt.savefig(path + str(rule)+scheme.replace(" ","") + 'ring'+str(n)+'Avg.pdf', bbox_inches='tight')
        if to_show:
            plt.show()
        plt.close()
        
        
def Graph_all_configurations(scheme, rules, df_magnetizacion, df_energia, camino, limites, fig_size, iteraciones = 1000, mostrar = False):
    n = 16
    path = camino + 'all\\'
    columnas = ["Esquema", "regla"] + [str(i) for i in range(iteraciones + 1)]
    df_magnetizacion = df_magnetizacion[columnas]
    df_energia = df_energia[columnas]
    if scheme == "block sequential":
        for r in rules:
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
            
            for B in df_mag_0.keys():
                fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
                ax1.set_ylim(limites[0])
                ax1.tick_params(axis='both', which='major', labelsize=5)
                ax1.set_ylabel("Density", fontsize=8)
                xx = [int(i) for i in df_mag_0[B].columns]
                for i in range(len(df_mag_0[B])):
                    p1 = ax1.plot(xx, df_mag_0[B].iloc[i], label = "Update Mode: "+str(i))                    
                plt.title("Rule "+str(r)+' '+scheme.title()+' ring '+str(n)+#"\n Average over Update Modes with "+str(B)+" blocks",
                           "\n with "+str(B)+" blocks",fontsize=8)
                plt.legend(loc="lower right", bbox_to_anchor=(.9, -.55), 
                            bbox_transform=fig.transFigure, ncol=3).set_visible(False)
                plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
                fig.set_size_inches(fig_size[0],fig_size[1])                
                plt.draw()
                plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'mag'+str(B)+'.pdf', bbox_inches='tight')
                if mostrar:
                    plt.show()
                plt.close()
                
                fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
                ax1.set_ylim(limites[1])
                ax1.tick_params(axis='both', which='major', labelsize=5)
                ax1.set_ylabel("Energy", fontsize=8)
                for i in range(len(df_mag_0[B])):
                    p2 = ax1.plot(xx, df_ene_0[B].iloc[i]/n, label = "Update Mode: "+str(i))#, linestyle='dotted')
                    
                plt.title("Rule "+str(r)+' '+scheme.title()+' ring '+str(n)+#"\n Average over Update Modes with "+str(B)+" blocks"
                            "\n with "+str(B)+" blocks"
                             , fontsize=8)
                plt.legend(loc="lower right", bbox_to_anchor=(.9, -.55), 
                            bbox_transform=fig.transFigure, ncol=3).set_visible(False)
                
                plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
                fig.set_size_inches(fig_size[0],fig_size[1])
                plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'ene'+str(B)+'.pdf', bbox_inches='tight')
                plt.draw()
                if mostrar:
                    plt.show()
                plt.close()
            
    elif scheme == 'sequential' or scheme == 'block parallel':
        for r in rules:
            df_mag_r = df_magnetizacion[df_magnetizacion["regla"] == str(r)+' '+scheme]
            df_mag_r = df_mag_r.reset_index().drop(columns = ["regla","Esquema","index"]).iloc[-32:].reset_index().drop(columns = ["index"])
            df_ene_r = df_energia[df_energia["regla"] == str(r)+' '+scheme]
            df_ene_r = df_ene_r.reset_index().drop(columns = ["regla","Esquema","index"]).iloc[-32:].reset_index().drop(columns = ["index"])
            
            fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
            ax1.set_ylim(limites[0]) 
            ax1.tick_params(axis='both', which='major', labelsize=5)
            ax1.set_ylabel("Density", fontsize=8)
            xx = [int(i) for i in df_ene_r.columns]
            for i in range(32):
                p1 = ax1.plot(xx, df_mag_r.iloc[i], label = "Update Mode: "+str(i))
            plt.title("Rule "+str(r)+' '+scheme.title()+' ring '+str(n), fontsize=8)#+"\n Average over different sample schemes.")
            plt.legend(loc="lower right", bbox_to_anchor=(1.1, -.5), 
                        bbox_transform=fig.transFigure, ncol=4).set_visible(False)
            plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
            fig.set_size_inches(fig_size[0],fig_size[1])
            fig.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'mag.pdf', bbox_inches='tight')
            if mostrar: 
                plt.show()
            plt.close()
            
            fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
            ax1.set_ylim(limites[1])
            ax1.tick_params(axis='both', which='major', labelsize=5)
            ax1.set_ylabel("Energy", fontsize=8)
            xx = [int(i) for i in df_ene_r.columns]
            for i in range(32):
                p2 = ax1.plot(xx, df_ene_r.iloc[i]/n, label = "Update Mode: "+str(i))#, linestyle='dotted')
            plt.title("Rule "+str(r)+' '+scheme.title()+' ring '+str(n), fontsize=8)#+"\n Average over different sample schemes.")
            plt.legend(loc="lower right", bbox_to_anchor=(1.1, -.5), 
                        bbox_transform=fig.transFigure, ncol=4).set_visible(False)
            plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
            fig.set_size_inches(fig_size[0],fig_size[1])
            plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'ene.pdf', bbox_inches='tight')
            if mostrar: 
                plt.show()
            plt.close()

    elif scheme == 'local clocks':
        P = [2,4,5]
        for r in rules:
            df_mag_r = df_magnetizacion[df_magnetizacion["regla"] == str(r)+' '+scheme]
            df_ene_r = df_energia[df_energia["regla"] == str(r)+' '+scheme]
            esquemas = [eval(elemento) for elemento in list(df_mag_r["Esquema"])]
            periodos = {str(esquema): max(esquema[1][0]) for esquema in esquemas}
            
            for periodo in P:
                aux = [str(esquema) for esquema in esquemas if periodos[str(esquema)] == periodo]
                df_mag_0 = df_mag_r[df_mag_r["Esquema"].isin(aux)].reset_index().drop(columns = ["regla","Esquema","index"])
                df_ene_0 = df_ene_r[df_ene_r["Esquema"].isin(aux)].reset_index().drop(columns = ["regla","Esquema","index"])
                
                fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
                ax1.set_ylim(limites[0])
                ax1.tick_params(axis='both', which='major', labelsize=5)
                ax1.set_ylabel("Density")#, fontsize=8)
                xx = [int(i) for i in df_mag_0.columns]
                for i in range(32):
                    ax1.plot(xx, df_mag_0.iloc[i], label = "Update Mode: "+str(i))
                plt.title("Rule "+str(r)+' '+scheme+' ring '+str(n)#+"\n Average over Update Modes 
                          +"\n with period "+str(periodo)
                           , fontsize=8)
                plt.legend(loc="lower right", bbox_to_anchor=(.9, -.55), 
                            bbox_transform=fig.transFigure, ncol=3).set_visible(False)
                plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
                fig.set_size_inches(fig_size[0],fig_size[1])
                plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'mag'+str(periodo)+'.pdf', bbox_inches='tight')
                if mostrar:
                    plt.show()
                plt.close()
                fig, ax1 = plt.subplots()#figsize=(8,5))#, layout='constrained')
                ax1.set_ylim(limites[1])
                ax1.tick_params(axis='both', which='major', labelsize=5)
                ax1.set_ylabel("Energy")#, fontsize=8)
                xx = [int(i) for i in df_mag_0.columns]
                for i in range(32):
                    ax1.plot(xx, df_ene_0.iloc[i]/n, label = "Update Mode: "+str(i))#, linestyle='dotted')
                plt.title("Rule "+str(r)+' '+scheme+' ring '+str(n)+#"\n Average over Update Modes with period "+str(periodo),
                           "\n with period "+str(periodo)
                            , fontsize=8)
                plt.legend(loc="lower right", bbox_to_anchor=(.9, -.55), 
                            bbox_transform=fig.transFigure, ncol=3).set_visible(False)
                plt.grid(color = 'gray', linestyle = '--', linewidth = 0.1)
                fig.set_size_inches(fig_size[0],fig_size[1])
                plt.savefig(path + str(r)+scheme.replace(" ","") + 'ring'+str(n)+'ene'+str(periodo)+'.pdf', bbox_inches='tight')
                if mostrar:
                    plt.show()
                plt.close()

def Graph_varianza(n, camino, mostrar = False):
    path = camino + 'varianza\\'
    esquemas = ["block sequential","local clocks","sequential", "block parallel"]
    reglas = [54, 90, 110, 150]
    df_ene = pd.read_csv("energia.csv", index_col=("Unnamed: 0"), sep = ",")
    df_mag = pd.read_csv("magnetizacion.csv", index_col=("Unnamed: 0"), sep = ",")
    df_ene = df_ene[df_ene["Tamano Muestra"] == 2**n].drop(columns = ["Tamano Muestra", "Largo Anillo"])
    df_mag = df_mag[df_mag["Tamano Muestra"] == 2**n].drop(columns = ["Tamano Muestra", "Largo Anillo"])

    varianzas_ene = {regla:{esquema:[] for esquema in esquemas} for regla in reglas}
    varianzas_mag = {regla:{esquema:[] for esquema in esquemas} for regla in reglas}

    esquema = 'block sequential'
    for r in reglas:
        df_ene_r = df_ene[df_ene["regla"] == str(r)+' '+esquema].drop(columns = ["regla"])
        df_mag_r = df_mag[df_mag["regla"] == str(r)+' '+esquema].drop(columns = ["regla"])
        
        mag_var_dic = {}
        mag_var_dic[3] = df_mag_r[~df_mag_r["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])
        varianzas_mag[r][esquema].append(np.var(np.mean(mag_var_dic[3][~mag_var_dic[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy(),axis=1)))
        mag_var_dic[3] = np.var(np.mean(mag_var_dic[3][~mag_var_dic[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy(), axis = 1))
        mag_var_dic[4] = df_mag_r[~df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
        varianzas_mag[r][esquema].append(np.var(np.mean(mag_var_dic[4][mag_var_dic[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy(),axis=1)))
        mag_var_dic[4] = np.var(np.mean(mag_var_dic[4][mag_var_dic[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy(), axis = 1))
        mag_var_dic[5] = np.var(np.mean(df_mag_r[df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"]).to_numpy(), axis = 1))
        varianzas_mag[r][esquema].append(np.var(np.mean(df_mag_r[df_mag_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"]).to_numpy(),axis=1)))

        ene_var_dic = {}
        ene_var_dic[3] = df_ene_r[~df_ene_r["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["index"])
        varianzas_ene[r][esquema].append(np.var(np.mean(ene_var_dic[3][~ene_var_dic[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy()/n)))
        ene_var_dic[3] = np.var(ene_var_dic[3][~ene_var_dic[3]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy()/n, axis = 1)
        ene_var_dic[4] = df_ene_r[~df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["index"])
        varianzas_ene[r][esquema].append(np.var(np.mean(ene_var_dic[4][ene_var_dic[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy()/n,axis=1)))
        ene_var_dic[4] = np.var(ene_var_dic[4][ene_var_dic[4]["Esquema"].str.contains("3: ")].reset_index().drop(columns = ["Esquema","index"])[-32:].to_numpy()/n, axis = 1)
        ene_var_dic[5] = np.var(df_ene_r[df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"]).to_numpy()/n, axis = 1)
        varianzas_ene[r][esquema].append(np.var(np.mean(df_ene_r[df_ene_r["Esquema"].str.contains("4: ")].reset_index().drop(columns = ["Esquema","index"]).to_numpy()/n,axis=1)))

        for i in [3,4,5]:
            plt.scatter(range(32),ene_var_dic[i])
            plt.hist(ene_var_dic[i])
            plt.title("Energy\nRule "+str(r)+' '+esquema+' ring '+str(n)+"\n Variance of different schemes with "+str(i)+" blocks")
            plt.savefig(path +'Varianza'+ str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'ene'+str(i)+'.pdf', bbox_inches='tight')
            if mostrar:
                plt.show()
            plt.close()
            # plt.scatter(range(32),mag_var_dic[i])
            plt.hist(mag_var_dic[i])
            plt.title("Density\nRule "+str(r)+' '+esquema+' ring '+str(n)+"\n Variance of different schemes with "+str(i)+" blocks")
            plt.savefig(path +'Varianza'+ str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'mag'+str(i)+'.pdf', bbox_inches='tight')
            if mostrar:
                plt.show()
            plt.close()

    esquema = 'sequential'
    for r in reglas:
        ene_var = np.var(df_ene[df_ene["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy()/n, axis = 1)
        varianzas_ene[r][esquema].append(np.var(np.mean(df_ene[df_ene["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy()/n,axis=1)))
        mag_var = np.var(df_mag[df_mag["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy(), axis = 1)
        varianzas_mag[r][esquema].append(np.var(np.mean(df_mag[df_mag["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy(),axis=1)))
        
        # plt.scatter(range(32),ene_var)
        plt.hist(ene_var)
        plt.title("Variance of Energy\nRule "+str(r)+' '+esquema+' ring '+str(n))
        plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'ene.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
        # plt.scatter(range(32),mag_var)
        plt.hist(mag_var)
        plt.title("Variance of Density\nRule "+str(r)+' '+esquema+' ring '+str(n))
        plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'mag.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
    esquema = 'block parallel'
    r = 110
    for r in reglas:
        ene_var = np.var(df_ene[df_ene["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy()/n, axis = 1)
        varianzas_ene[r][esquema].append(np.var(np.mean(df_ene[df_ene["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy()/n,axis = 1)))
        mag_var = np.var(df_mag[df_mag["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy(), axis = 1)
        varianzas_mag[r][esquema].append(np.var(np.mean(df_mag[df_mag["regla"] == str(r)+' '+esquema].reset_index().drop(columns = ["regla","Esquema","index"])[-32:].to_numpy(),axis=1)))
        
        # plt.scatter(range(32),ene_var)
        plt.hist(ene_var)
        plt.title("Variance of Energy\nRule "+str(r)+' '+esquema+' ring '+str(n))
        plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'ene.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
        # plt.scatter(range(32),mag_var)
        plt.hist(mag_var)
        plt.title("Variance of Density\nRule "+str(r)+' '+esquema+' ring '+str(n))
        plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'mag.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
    esquema = 'local clocks'
    P = [2,4,5]
    for r in reglas:
        df_mag_r = df_mag[df_mag["regla"] == str(r)+' '+esquema]
        df_ene_r = df_ene[df_ene["regla"] == str(r)+' '+esquema]
        
        esquemass = [eval(elemento) for elemento in list(df_ene_r["Esquema"])]
        periodos = {str(esquema): max(esquema[1][0]) for esquema in esquemass}
        for periodo in P:
            aux = [str(esquema) for esquema in esquemass if periodos[str(esquema)] == periodo]
            mag_var = np.var(df_mag_r[df_mag_r["Esquema"].isin(aux)].reset_index().drop(columns = ["Esquema","regla","index"])[-32:].to_numpy(),axis = 1)
            varianzas_mag[r][esquema].append(np.var(np.mean(df_mag_r[df_mag_r["Esquema"].isin(aux)].reset_index().drop(columns = ["Esquema","regla","index"])[-32:].to_numpy(),axis = 1)))
            fig = plt.figure()
            # plt.scatter(range(32),mag_var)
            plt.hist(mag_var)
            plt.title("Density\nRule "+str(r)+' '+esquema+' ring '+str(n)+"\n Variance of different schemes with period "+str(periodo))
            plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") +'Ring'+str(n)+'mag'+str(periodo)+'.pdf', bbox_inches='tight')
            if mostrar:
                plt.show()
            plt.close()
            ene_var = np.var(df_ene_r[df_ene_r["Esquema"].isin(aux)].reset_index().drop(columns = ["Esquema","regla","index"])[-32:].to_numpy()/n,axis = 1)
            varianzas_ene[r][esquema].append(np.var(np.mean(df_ene_r[df_ene_r["Esquema"].isin(aux)].reset_index().drop(columns = ["Esquema","regla","index"])[-32:].to_numpy()/n,axis=1)))
            fig = plt.figure()
            # plt.scatter(range(32),ene_var)    
            plt.hist(ene_var)
            plt.title("Energy\nRule "+str(r)+' '+esquema+' ring '+str(n)+"\n Variance of different schemes with period "+str(periodo))
            plt.savefig(path +'Varianza' + str(r)+esquema.replace(" ","") + 'Ring'+str(n)+'ene'+str(periodo)+'.pdf', bbox_inches='tight')
            if mostrar:
                plt.show()
            plt.close()
            
    etiquetas = ["sequential", "block sequential 3","block sequential 4",
                 "block sequential 5","block parallel",
                 "local clocks 2", "local clocks 4", "local clocks 5"]
    for r in reglas:
        yy = []
        yy.append(varianzas_ene[r]["sequential"][0])
        for i in range(3):
            yy.append(varianzas_ene[r]["block sequential"][i])
        
        yy.append(varianzas_ene[r]["block parallel"][0])
        for i in range(3):
            yy.append(varianzas_ene[r]["local clocks"][i])    
        plt.scatter(range(len(yy)),yy)
        plt.xticks(range(len(yy)), etiquetas,rotation=80) 
        plt.grid(color='k', linestyle='-', linewidth=.2)
        plt.title("Variance of Energy\nRule "+str(r)+' ring '+str(n))
        plt.savefig(path + str(r)+'varianzaRing'+str(n)+'ene.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
        yy = []
        yy.append(varianzas_mag[r]["sequential"][0])
        for i in range(3):
            yy.append(varianzas_mag[r]["block sequential"][i])
        
        yy.append(varianzas_mag[r]["block parallel"][0])
        for i in range(3):
            yy.append(varianzas_mag[r]["local clocks"][i])    
        plt.scatter(range(len(yy)),yy)
        plt.xticks(range(len(yy)), etiquetas,rotation=80) 
        plt.grid(color='k', linestyle='-', linewidth=.2)
        plt.title("Variance of Density\nRule "+str(r)+' ring '+str(n))
        plt.savefig(path + str(r)+'varianzaRing'+str(n)+'mag.pdf', bbox_inches='tight')
        if mostrar:
            plt.show()
        plt.close()
    
def obtenerPyDelta(df, largo_anillo):
    # if "Tamano Muestra" in df.columns:
        # df = df[df["Tamano Muestra"] == 2**largo_anillo].drop(columns = ["Tamano Muestra", "Largo Anillo"])
    df = df[df["regla"].str.contains("local clocks")].reset_index().drop(columns = ["index"])
    lista = [eval(elemento) for elemento in list(df["Esquema"])]
    df = df.drop(columns = ["Esquema"])
    
    LC = {i:elemento[1] for i,elemento in enumerate(lista)}
    
    aux = []
    for k in LC.keys():
    # k = 2
        reloj = []
        fase = []
        try:
            for i in range(largo_anillo):
                bandera = False        
                for llave in LC[k].keys():
                    if len(LC[k]) == 2:
                        if i in LC[k][0]:
                            fase.append(0)
                            if i in LC[k][1]:
                                reloj.append(1)
                            else:
                                reloj.append(2)
                        else:
                            fase.append(1)
                            reloj.append(2)
                        break
                    else:
                        if bandera and i in LC[k][llave]:
                            reloj.append(llave-fase[i])
                            break
                        if not bandera and i in LC[k][llave]:
                            fase.append(llave)
                            bandera = True
            aux.append(["local clocks",[reloj,fase]])
        except AttributeError:
            aux.append(["local clocks",LC[k]])
                
    df["Esquema"] = aux
    
    cols = df.columns.tolist()
    cols = cols[-1:]+cols[:-1]
    df = df[cols]
    df = df.reset_index().drop(columns = ["index"])
    return df