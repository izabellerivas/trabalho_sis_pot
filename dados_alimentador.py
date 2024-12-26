from itertools import combinations
import numpy as np

rcalcJ = True
e = 0.00001
max_ite =30
n = 5
slack = 1 
PV = [4]
PQ = [2,3,5]
Sb = 100.0 #MVA
Vb = 13.8 #kV

P = [0.0]*n #vetor das potências ativas
Q = [0.0]*n #vetor das potências reativas

aux = 1
V = [complex(aux)]*n #vetor de tensão 
V[slack-1] = 1.02 + 0j

pares_de_barras = list(combinations(range(1, n+1), 2))
#Dados das barras PV
if PV != None:
    aux_p_PV = [(80.0-25.0)] 
    aux_v_PV = [1.0] 
    for i,id in enumerate(PV):
        P[id-1] = aux_p_PV[i]/Sb
        V[id-1] = complex(aux_v_PV[i])

#Dados das barras PQ
if PQ != None:
    aux_p_PQ = [-45,-25,-55]
    aux_q_PQ = [-15,-15,-20]
    for i,id in enumerate(PQ):
        P[id-1] = aux_p_PQ[i]/Sb
        Q[id-1] = aux_q_PQ[i]/Sb

admitancias={(1,1):0,
            (1,2):1/(0.02+0.2j),
            (1,3):0,
            (1,4):0,
            (1,5):1/(0.015+0.04j),
            (2,2):0,
            (2,3):1/(0.01+0.025j),
            (2,4):0,
            (2,5):0,
            (3,3):0,
            (3,4):1/(0.02+0.4j),
            (3,5):1/(0.02+0.05j),
            (4,4):0,
            (4,5):1/(0.015+0.04j),
            (5,5):0}

Y = np.zeros((n,n),dtype=complex)
for i in range(n):
    for j in range(n):
        try:
            Y[i,i] += admitancias[(i+1,j+1)]
            if i != j: Y[i,j] = - admitancias[(i+1,j+1)] 
        except:
            Y[i,i] += admitancias[(j+1,i+1)]
            Y[i,j] = - admitancias[(j+1,i+1)]

np.set_printoptions(linewidth=200)
print("Matriz de Admitância (em pu)")
print(Y)
print("\n")