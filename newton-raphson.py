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

def calc_dPdV(Y,V,Fase,id_i,id_j):
    output = 0
    i = id_i-1
    j = id_j-1
    if i == j:
        for k in range(len(V)):
            output += abs(V[k]*Y[i,k])*math.cos(cmath.phase(Y[i,k])-Fase[k]+Fase[k]) if i!= k else 2*abs(V[i]*Y[i,i])*math.cos(cmath.phase(Y[i,i]))
        return output
    else:
        return abs(V[i]*Y[i,j])*math.cos(cmath.phase(Y[i,j])-Fase[j]+Fase[j])

def calc_dQddelta(Y,V,Fase,id_i,id_j):
    output = 0
    i = id_i-1
    j = id_j-1
    if i == j:
        for k in range(len(V)):
            output += abs(V[i]*V[k]*Y[i,k])*math.cos(cmath.phase(Y[i,k])-Fase[i]+Fase[k]) if i!= k else 0
        return output
    else:
        return -abs(V[i]*V[j]*Y[i,j])*math.cos(cmath.phase(Y[i,j])-Fase[i]+Fase[j])

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

np.set_printoptions(linewidth=200)

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

Pinj_PQ = [pot for i,pot in enumerate(dados.P) if (i+1) in dados.PQ]
Pinj_PV = [pot for i,pot in enumerate(dados.P) if (i+1) in dados.PV]
Qinj_PQ = [potq for i,potq in enumerate(dados.Q) if (i+1) in dados.PQ]

inj = Pinj_PQ + Pinj_PV + Qinj_PQ
V = [abs(num) for num in dados.V]
Fase = [cmath.phase(num) for num in dados.V]
V_ite = [V]
Fase_ite = [Fase]

P = dados.P.copy()
Q = dados.Q.copy()
for id in dados.PQ:
    P[id-1] = calc_P(dados.Y,V,Fase,id)
    Q[id-1] = calc_Q(dados.Y,V,Fase,id)
for id in dados.PV:
    P[id-1] = calc_P(dados.Y,V,Fase,id)

Q_ite = [Q]
P_ite = [P]

erro = 1
cont = 0
PQPV = dados.PQ + dados.PV
PQ_PV_PQ = dados.PQ + dados.PV + dados.PQ 

while (cont < dados.max_ite):
    novoV = V_ite[-1].copy()
    novaFase = Fase_ite[-1].copy()
    novoQ = Q_ite[-1].copy()
    novoP = P_ite[-1].copy()

    for id in PQPV:
        novoP[id-1] = calc_P(dados.Y,novoV,novaFase,id)
        if id in dados.PQ:
            novoQ[id-1] = calc_Q(dados.Y,novoV,novaFase,id)

    novoPinj_PQ = [pot for i,pot in enumerate(novoP) if (i+1) in dados.PQ]
    novoPinj_PV = [pot for i,pot in enumerate(novoP) if (i+1) in dados.PV]
    novoQinj_PQ = [potq for i,potq in enumerate(novoQ) if (i+1) in dados.PQ]
    novoInj = novoPinj_PQ+novoPinj_PV+novoQinj_PQ
    # print("novoP: ",novoP)
    # print("novoQ: ",novoQ)
    # print("novoInj: ",novoInj)
    # print("inj certo: ",inj)
    D = np.array([base - novo for base,novo in zip(inj,novoInj)]) 

    erro = max(np.vectorize(abs)(D)) #precisa aplicar o abs em cada termo antes
    if (erro<dados.e):
        V_ite += [novoV]
        P_ite += [novoP]
        Q_ite += [novoQ]
        Fase_ite += [novaFase] 
        break
    cont += 1
    print(f"\nIteração {cont}")

    if dados.rcalcJ or (cont == 0):
        for i in range(nJ):
            for j in range(nJ):
                id_i = PQ_PV_PQ[i]
                id_j = PQ_PV_PQ[j]

                if (i < nPQ + nPV) and (j < nPQ + nPV):
                    J[i,j] = calc_dPddelta(dados.Y,novoV,novaFase,id_i,id_j)
                elif (i < nPQ + nPV) and (j >= nPQ + nPV):
                    J[i,j] = calc_dPdV(dados.Y,novoV,novaFase,id_i,id_j)
                elif (i >= nPQ + nPV) and (j < nPQ + nPV):
                    J[i,j] = calc_dQddelta(dados.Y,novoV,novaFase,id_i,id_j)
                elif (i >= nPQ + nPV) and (j >= nPQ + nPV):
                    J[i,j] = calc_dQdV(dados.Y,novoV,novaFase,id_i,id_j)

        # print("Jacobiano:")
        # print(J)
        invJ = np.linalg.inv(J)
        # print("Jacobiano^-1:")
        # print(invJ)

    #cálculo do erro da potencia
    # R = [J]^-1*D
    
    print("D:",D)
    R = invJ @ D
    print("R: ",R)
    auxPQ = 0
    for i,id in enumerate(PQPV):
        novaFase[id-1] += R[i]
        novoP[id-1] = calc_P(dados.Y,novoV,novaFase,id)
        if id in dados.PQ:
            novoV[id-1] += R[nPQ+nPV+auxPQ]
            novoQ[id-1] = calc_Q(dados.Y,novoV,novaFase,id)
            auxPQ += 1
            
    # print("NovoV: ",novoV)
    # print("NovaFase: ",novaFase)
    # print("novoP: ",novoP)
    # print("novoQ: ",novoQ)

    erro = max(np.vectorize(abs)(D)) #precisa aplicar o abs em cada termo antes
    V_ite += [novoV]
    P_ite += [novoP]
    Q_ite += [novoQ]
    Fase_ite += [novaFase]
    
    print([f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(novoV,novaFase)],"\terro: ", erro)

print("\nResultados de tensão por iteração")
for i,lista in enumerate(V_ite):
    print(f"{i}: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(lista,Fase_ite[i])])

print("\nResultados de potência por iteração")
for i,lista in enumerate(P_ite):
    print(f"{i}: ",[f"{p:.4f}+ j({q:.4f}) pu" for p,q in zip(lista,Q_ite[i])])

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
print(f"Número de iterações: {cont}")
print("V: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in zip(novoV,novaFase)])
print("S: ",[f"{p:.4f}+ j({q:.4f}) pu" for p,q in zip(novoP,novoQ)])
print("S: ",[f"{p*dados.Sb:.4f}+ j({q*dados.Sb:.4f}) MVA" for p,q in zip(novoP,novoQ)])
