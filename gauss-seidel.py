'''Gauss-Seidel aplicado para cálculo do fluxo de potência'''

from itertools import combinations
import numpy as np
import cmath
import math
import dados as dados

np.set_printoptions(linewidth=200)

#Método Iterativo - Gauss-Seidel
V_ite = [dados.V]
Q_ite = [dados.Q]
erro = 1
cont = 0
print("it: V | Q | erro")
while (erro > dados.e) and (cont < dados.max_ite):
    cont += 1
    novoQ = Q_ite[-1].copy()
    novoV = V_ite[-1].copy()
    if dados.PV != None:
        for id in dados.PV:
            i = id-1
            tempQ = 0
            for j in range(dados.n):
                tempQ -= abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(dados.Y[i,j])-cmath.phase(novoV[i])+cmath.phase(novoV[j]))
            novoQ[i] = tempQ
        for i in range(dados.n):
            if i != dados.slack-1:
                aux = (complex(dados.P[i],-novoQ[i])/novoV[i].conjugate())
                for j in range(dados.n):
                    aux -= dados.Y[i,j]*novoV[j] if i!=j else 0
                temp = aux/dados.Y[i,i]
                if (i+1) in dados.PV:
                    nfase = cmath.phase(complex(math.sqrt(abs(dados.V[i])*abs(dados.V[i]) - temp.imag*temp.imag),temp.imag))
                    novoV[i] = cmath.rect(abs(dados.V[i]),nfase)
                else: 
                    novoV[i] = temp
        aux_erro = max([abs(abs(a)- abs(b)) for a, b in zip(novoV, V_ite[-1])])
        erro = max(aux_erro,max([abs(abs(a) - abs(b)) for a, b in zip(novoQ, Q_ite[-1])])) if cont !=1 else 1
        print(f"{cont}: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in list(map(cmath.polar,novoV))],[f"{valor:.4f}" for valor in novoQ],"\terro: ", erro)
        V_ite += [novoV]
        Q_ite += [novoQ]
    else:
        for i in range(dados.n):
            if i != dados.slack-1:
                aux = (complex(dados.P[i],-novoQ[i])/novoV[i].conjugate())
                for j in range(dados.n):
                    aux -= dados.Y[i,j]*novoV[j] if i!=j else 0
                novoV[i] = aux/dados.Y[i,i]
        erro = max([abs(a)- abs(b) for a, b in zip(novoV, V_ite[-1])])
        print(f"{cont}: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in list(map(cmath.polar,novoV))], "\terro: ", erro)
        V_ite += [novoV]

print("\n")
#Cálculo da potência
P_l = np.zeros((dados.n,dados.n))
Q_l = np.zeros((dados.n,dados.n))
P_final = [0]*dados.n
Q_final = [0]*dados.n
print("Potência injetada na barra")
for i in range(dados.n):
    for j in range(dados.n):
        if (i+1) in dados.PQ:
            P_final[i] = dados.P[i]
            Q_final[i] = dados.Q[i]
        elif (i+1) in dados.PV:
            P_final[i] = dados.P[i]
            Q_final[i] = novoQ[i]
        else:
            P_final[i] += abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.cos(cmath.phase(dados.Y[i,j])-cmath.phase(novoV[i])+cmath.phase(novoV[j]))
            Q_final[i] -= abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(dados.Y[i,j])-cmath.phase(novoV[i])+cmath.phase(novoV[j]))
            
        P_l[i,j] = 0 if i==j else -abs(dados.Y[i,j]*novoV[i]*novoV[i])*math.cos(cmath.phase(dados.Y[i,j]))+abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.cos(cmath.phase(dados.Y[i,j])-cmath.phase(novoV[i])+cmath.phase(novoV[j]))
        Q_l[i,j] = 0 if i==j else abs(dados.Y[i,j]*novoV[i]*novoV[i])*math.sin(cmath.phase(dados.Y[i,j]))-abs(dados.Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(dados.Y[i,j])-cmath.phase(novoV[i])+cmath.phase(novoV[j]))

    print(f"Barra {i+1}: {P_final[i]:.4f} + j({Q_final[i]:.4f}) pu \t\t {P_final[i]*dados.Sb:.4f} + j({Q_final[i]*dados.Sb:.4f}) MVA ")

print("\nMatriz de potência ativa nas linhas (em pu)")
print(P_l)
print("\nMatriz de potência reativa nas linhas (em pu)")
print(Q_l)

print("\n")

perdasP = {}
perdasQ = {}
for barra1, barra2 in dados.pares_de_barras:
    perdasP[(barra1,barra2)] = P_l[barra1-1,barra2-1]+ P_l[barra2-1,barra1-1]
    perdasQ[(barra1,barra2)] = Q_l[barra1-1,barra2-1]+ Q_l[barra2-1,barra1-1]
    print(f"Perdas na linha {barra1}-{barra2}: {perdasP[(barra1,barra2)]:.4f} + j({perdasQ[(barra1,barra2)]:.4f}) pu \t\t {perdasP[(barra1,barra2)]*dados.Sb:.4f} + j({perdasQ[(barra1,barra2)]*dados.Sb:.4f}) MVA")

perdastot_P = np.sum(P_l)
perdastot_Q = np.sum(Q_l)

print(f"Número de iterações: {cont}")
print(f"\nPerdas totais: {perdastot_P:.4f} + j({perdastot_Q:.4f}) pu \t\t {perdastot_P*dados.Sb:.4f} + j({perdastot_Q*dados.Sb:.4f}) MVA")

print("V: ",[f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in list(map(cmath.polar,novoV))])
print("S: ",[f"{p:.4f}+ j({q:.4f}) pu" for p,q in zip(P_final,Q_final)])
print("S: ",[f"{p*dados.Sb:.4f}+ j({q*dados.Sb:.4f}) MVA" for p,q in zip(P_final,Q_final)])