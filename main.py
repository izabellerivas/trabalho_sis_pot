from itertools import combinations
import numpy as np
import cmath
import math

flag = True #dados adicionados diretamente no script

if not flag:
    #coleta de dados iniciais
    print("Caso algum valor não exista, escreva -")
    e = float(input("Erro para convergência: "))
    n = int(input("Número de barras: "))
    slack = int(input("Qual o id da barra slack? (Iniciando no 1) "))

    aux = input("Quais os ids das barras PV? (ex: 2 3 4) ")
    PV = None if aux == '-' else [int(num) for num in aux.split()]

    aux = input("Quais os ids das barras PQ? (ex: 2 3 4) ")
    PQ = None if aux == '-' else [int(num) for num in aux.split()]

    Sb = float(input("Qual a potência base? (em MVA) "))
    Vb = float(input("Qual a tensão base? (em kV) "))

    P = [0.0]*n #vetor das potências ativas
    Q = [0.0]*n #vetor das potências reativas

    aux = input("Qual o chute inicial de tensão das barras? ")
    V = [float(aux)]*n #vetor de tensão 
    V[slack-1] = complex((input("Qual a tensão da barra slack? é assumido fase igual a 0 (em pu) ")))


    print("\nColeta de dados de potência, caso haja carga e geração conectadas insira apenas a líquida\n")
    #Dados das barras PV
    if PV != None:
        aux_p_PV = input("Potência das barras PV? (insira em MW como lista na mesma ordem inserida para os ids) ").split()
        aux_v_PV = input("Tensão das barras PV? (insira em pu como lista na mesma ordem inserida para os ids) ").split()
        for i,id in enumerate(PV):
            P[id-1] = float(aux_p_PV[i])/Sb
            V[id-1] = float(aux_v_PV[i])
    #Dados das barras PQ
    if PQ != None:
        aux_p_PQ = input("Potência ativa das barras PQ? (insira em MW como lista na mesma ordem inserida para os ids) ").split()
        aux_q_PQ = input("Potência reativa das barras PQ? (insira em MVAr como lista na mesma ordem inserida para os ids) ").split()
        for i,id in enumerate(PQ):
            P[id-1] = float(aux_p_PQ[i])/Sb
            Q[id-1] = float(aux_q_PQ[i])/Sb

    #CRIAÇÃO DA MATRIZ DE ADMITÂNCIA
    print("\nColeta da reatância das linhas\n")
    print("Caso uma barra não esteja conectada na outra coloque -")
    Zb = Vb*Vb/Sb
    admitancias = {}
    pares_de_barras = list(combinations(range(1, n+1), 2))

    for barra1, barra2 in pares_de_barras:
        aux = input(f"Informe a reatância entre a barra {barra1} e a barra {barra2} em pu: (ex: 1+3j) ")
        admitancia = 0 if aux == '-' else 1/complex(aux)
        admitancias[(barra1, barra2)] = admitancia
    for i in range(n):
        aux = input(f"Informe a reatância da barra {i+1} em pu: ")
        admitancia = 0 if aux == '-' else 1/complex(aux)
        admitancias[(i+1,i+1)] = admitancia
else:
    e = 0.01
    n = 4 
    slack = 1 
    PV = None
    PQ = [2,3,4]
    Sb = float(100.0) #MVA
    Vb = float(13.8) #kV

    delta = [0]*n #ângulo das tensões
    P = [0.0]*n #vetor das potências ativas
    Q = [0.0]*n #vetor das potências reativas

    aux = 1
    V = [complex(aux)]*n #vetor de tensão 
    V[slack-1] = 1.05 + 0j

    #Dados das barras PV
    if PV != None:
        aux_p_PV = ['50','100','150'] 
        aux_v_PV = ['1','1','1'] 
        for i,id in enumerate(PV):
            P[id-1] = float(aux_p_PV[i])/Sb
            V[id-1] = float(aux_v_PV[i])

    #Dados das barras PQ
    if PQ != None:
        aux_p_PQ = ['50','100','150'] 
        aux_q_PQ = ['30','45','-105']
        for i,id in enumerate(PQ):
            P[id-1] = float(aux_p_PQ[i])/Sb
            Q[id-1] = float(aux_q_PQ[i])/Sb

    admitancias={(1,1):1/(1.25j),
                (1,2):0,
                (1,3):1/(0.25j),
                (1,4):1/(0.2j),
                (2,2):1/(1.25j),
                (2,3):1/(0.4j),
                (2,4):1/(0.2j),
                (3,3):1/(1.25j),
                (3,4):1/(0.125j),
                (4,4):0}

Y = np.zeros((n,n),dtype=complex)
for i in range(n):
    for j in range(n):
        try:
            Y[i,i] += admitancias[(i+1,j+1)]
            if i != j: Y[i,j] = - admitancias[(i+1,j+1)] 
        except:
            Y[i,i] += admitancias[(j+1,i+1)]
            Y[i,j] = - admitancias[(j+1,i+1)]
print("Matriz de Admitância (em pu)")
print(Y)
print("\n")
#Método Iterativo - Gauss-Seidel
V_ite = [V]
Q_ite = [Q]
erro = 1
while erro > e:
    novoQ = Q_ite[-1].copy()
    novoV = V_ite[-1].copy()
    if PV != None:
        for id in PQ:
            novoQ[id-1] = 0
    else:
        for i in range(n):
            if i != slack-1:
                aux = (complex(P[i],-novoQ[i])/novoV[i].conjugate())
                for j in range(n):
                    aux -= Y[i,j]*novoV[j] if i!=j else 0
                novoV[i] = aux/Y[i,i]
        # print([f"{abs(x):.5f}" for x in novoV])
        erro = max([abs(a - b) for a, b in zip(novoV, V_ite[-1])])
        print([f"{modulo:.4f}<{fase*180/cmath.pi:.2f}°" for modulo,fase in list(map(cmath.polar,novoV))], "\terro: ", erro)
        V_ite += [novoV]
        # print(f"Erro: {erro}")
print("\n")
#Cálculo da potência
P_l = np.zeros((n,n))
Q_l = np.zeros((n,n))
P_final = P.copy()
Q_final = Q.copy()
for i in range(n):
    for j in range(n):
        P_final[i] += abs(Y[i,j]*novoV[i]*novoV[j])*math.cos(cmath.phase(Y[i,j])-cmath.phase(novoV[i]+cmath.phase(novoV[j])))
        Q_final[i] -= abs(Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(Y[i,j])-cmath.phase(novoV[i]+cmath.phase(novoV[j])))

        P_l[i,j] = 0 if i==j else abs(Y[i,j]*novoV[i]*novoV[j])*math.cos(cmath.phase(Y[i,j])-cmath.phase(novoV[i]+cmath.phase(novoV[j])))
        Q_l[i,j] = 0 if i==j else -abs(Y[i,j]*novoV[i]*novoV[j])*math.sin(cmath.phase(Y[i,j])-cmath.phase(novoV[i]+cmath.phase(novoV[j])))

    print(f"Barra {i+1}: {P_final[i]:.4f} + j({Q_final[i]:.4f}) pu \t\t {P_final[i]*Sb:.4f} + j({Q_final[i]*Sb:.4f}) MVA ")

print("\nMatriz de potência ativa nas linhas (em pu)")
print(P_l)
print("\nMatriz de potência reativa nas linhas (em pu)")
print(Q_l)

perdas_P = np.sum(P_l)
perdas_Q = np.sum(Q_l)
print(f"\nPerdas: {perdas_P} + j({perdas_Q}) pu \t\t {perdas_P*Sb} + j({perdas_Q*Sb}) MVA")