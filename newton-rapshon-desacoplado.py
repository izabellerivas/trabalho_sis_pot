from itertools import combinations
import numpy as np
import math
import cmath
import dados_alimentador as dados

#Calcula a potência ativa na barra id (lembrando que o id começa no 1)
def calc_P(Y,V,Fase,id):
    P = 0
    i = id - 1
    for j in range(len(V)):
        P += abs(V[i]*V[j]*Y[i,j])*math.cos(cmath.phase(Y[i,j])-Fase[i]+Fase[j])
    return P

#Calcula a potência reativa na barra id (lembrando que o id começa no 1)
def calc_Q(Y,V,Fase,id):
    Q = 0
    i = id - 1
    for j in range(len(V)):
        Q += -abs(V[i]*V[j]*Y[i,j])*math.sin(cmath.phase(Y[i,j])-Fase[i]+Fase[j])
    return Q

def calc_dPddelta(Y,V,Fase,id_i,id_j):
    output = 0
    i = id_i-1
    j = id_j-1
    if i == j:
        for k in range(len(V)):
            output += abs(V[i]*V[k]*Y[i,k])*math.sin(cmath.phase(Y[i,k])-Fase[k]+Fase[k]) if i!= k else 0
        return output
    else:
        return -abs(V[i]*V[j]*Y[i,j])*math.sin(cmath.phase(Y[i,j])-Fase[j]+Fase[j])


def calc_dQdV(Y,V,Fase,id_i,id_j):
    output = 0
    i = id_i-1
    j = id_j-1
    if i == j:
        for k in range(len(V)):
            output += -abs(V[k]*Y[i,k])*math.sin(cmath.phase(Y[i,k])-Fase[i]+Fase[k]) if i!= k else -2*abs(V[i]*Y[i,i])*math.sin(cmath.phase(Y[i,i]))
        return output
    else:
        return -abs(V[i]*Y[i,j])*math.sin(cmath.phase(Y[i,j])-Fase[i]+Fase[j])

print("Matriz de Admitância (em pu)")
print(dados.Y)
print("\n")

#NEWTON-RAPHSON
#montando matriz jacobiana
#barra PQ -> linha delta P e delta Q
#barra PV -> linha delta Q
#J -> número de barras PQ*2 + número de barras PV x número de barras PQ*2
nPQ = len(dados.PQ if dados.PQ is not None else [])
nPV = len(dados.PV if dados.PV is not None else [])
nJ = 2*nPQ + nPV
J = np.zeros((nJ,nJ))
R = np.zeros(nJ) #vetor de valores do delta fase e módulos da tensão
D = np.zeros(nJ) # vetor de valor do erro de potência

Pinj = [pot for i,pot in enumerate(dados.P) if i != (dados.slack-1)]
Qinj = [potq for i,potq in enumerate(dados.Q) if (i+1) in dados.PQ]

inj = Pinj + Qinj

V_ite = [[abs(num) for num in dados.V]]
Fase_ite = [[cmath.phase(num) for num in dados.V]]
Q_ite = [dados.Q]
P_ite = [dados.P]

erro = 1
cont = 0
PQPV = dados.PQ + dados.PV
PQ_PV_PQ = dados.PQ + dados.PV + dados.PQ 

while (erro > dados.e) | (cont > dados.max_ite):

    novoV = V_ite[-1].copy()
    novaFase = Fase_ite[-1].copy()
    novoQ = Q_ite[-1].copy()
    novoP = P_ite[-1].copy()

    if dados.rcalcJ or (cont == 0):
        for i in range(nJ):
            for j in range(nJ):
                id_i = PQ_PV_PQ[i]
                id_j = PQ_PV_PQ[j]

                if (i < nPQ + nPV) and (j < nPQ + nPV):
                    J[i,j] = calc_dPddelta(dados.Y,novoV,novaFase,id_i,id_j)
                elif (i < nPQ + nPV) and (j >= nPQ + nPV):
                    J[i,j] = 0
                elif (i >= nPQ + nPV) and (j < nPQ + nPV):
                    J[i,j] = 0
                elif (i >= nPQ + nPV) and (j >= nPQ + nPV):
                    J[i,j] = calc_dQdV(dados.Y,novoV,novaFase,id_i,id_j)

        print("Jacobiano:")
        print(J)
        invJ = np.linalg.inv(J)
        print("Jacobiano^-1:")
        print(invJ)

    #cálculo do erro da potencia
    # R = [J]^-1*D

    novoPinj = [pot for i,pot in enumerate(novoP) if i != (dados.slack-1)]
    novoQinj = [potq for i,potq in enumerate(novoQ) if (i+1) in dados.PQ]
    novoD = novoPinj+novoQinj
    print("novoP: ",novoP)
    print("novoQ: ",novoQ)
    print("novoD: ",novoD)
    D = np.array([base - novo for base,novo in zip(inj,novoD)]) if cont != 0 else np.array(novoD)
    print("D:",D)
    R = invJ @ D

    for i,id in enumerate(PQPV):
        novaFase[id-1] += R[i]
        novoP[id-1] = calc_P(dados.Y,novoV,novaFase,id)
        if id in dados.PQ:
            novoV[id-1] += R[nPQ+nPV+i]
            novoQ[id-1] = calc_Q(dados.Y,novoV,novaFase,id)
            

    erro = max(np.vectorize(abs)(D)) #precisa aplicar o abs em cada termo antes
    V_ite += [novoV]
    Q_ite += [novoQ]
    P_ite += [novoP]
    Q_ite += [novoQ]
    Fase_ite += [novaFase]
    cont += 1
    
    print([f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(novoV,novaFase)],[f"{valor:.4f}" for valor in novoQ],"\terro: ", erro)

print("\nResultados de tensão por iteração")
for i,lista in enumerate(V_ite):
    print(f"{i}: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(lista,Fase_ite[i])])

print("\nResultados de potência por iteração")
for i,lista in enumerate(P_ite):
    print([f"{p:.4f}+ j({q:.4f}) pu" for p,q in zip(lista,Q_ite[i])])

#Valores finais
#barra slack
novoP[dados.slack-1] = calc_P(dados.Y,novoV,novaFase,dados.slack)
novoQ[dados.slack-1] = calc_Q(dados.Y,novoV,novaFase,dados.slack)

#barra PV
if dados.PV != None:
    for id in dados.PV:
        novoQ[id-1] = calc_Q(dados.Y,novoV,novaFase,id)

#Potência nas linhas

P_l = np.zeros((dados.n,dados.n))
Q_l = np.zeros((dados.n,dados.n))

for i in range(dados.n):
    for j in range(dados.n):
        P_l[i,j] = 0 if i==j else -abs(dados.Y[i,j]*novoV[i]*novoV[i])*math.cos(cmath.phase(dados.Y[i,j]))+abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.cos(cmath.phase(dados.Y[i,j])-novaFase[i]+novaFase[j])
        Q_l[i,j] = 0 if i==j else abs(dados.Y[i,j]*novoV[i]*novoV[i])*math.sin(cmath.phase(dados.Y[i,j]))-abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(dados.Y[i,j])-novaFase[i]+novaFase[j])

print("\nMatriz de potência ativa nas linhas (em pu)")
print(P_l)
print("\nMatriz de potência reativa nas linhas (em pu)")
print(Q_l)

print("\n")

perdasP = {}
perdasQ = {}
print("\nResultados Finais")
for barra1, barra2 in dados.pares_de_barras:
    perdasP[(barra1,barra2)] = P_l[barra1-1,barra2-1]+ P_l[barra2-1,barra1-1]
    perdasQ[(barra1,barra2)] = Q_l[barra1-1,barra2-1]+ Q_l[barra2-1,barra1-1]
    print(f"Perdas na linha {barra1}-{barra2}: {perdasP[(barra1,barra2)]:.4f} + j({perdasQ[(barra1,barra2)]:.4f}) pu \t\t {perdasP[(barra1,barra2)]*dados.Sb:.4f} + j({perdasQ[(barra1,barra2)]*dados.Sb:.4f}) MVA")

perdastot_P = np.sum(P_l)
perdastot_Q = np.sum(Q_l)

print(f"\nPerdas totais: {perdastot_P:.4f} + j({perdastot_Q:.4f}) pu \t\t {perdastot_P*dados.Sb:.4f} + j({perdastot_Q*dados.Sb:.4f}) MVA")


#printando os resultados


print("V: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(novoV,novaFase)])
print("S: ",[f"{p:.4f}+ j({q:.4f}) pu" for p,q in zip(novoP,novoQ)])
print("S: ",[f"{p*dados.Sb:.4f}+ j({q*dados.Sb:.4f}) MVA" for p,q in zip(novoP,novoQ)])
AQUI = 1