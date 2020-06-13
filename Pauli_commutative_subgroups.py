import numpy as np

def tobinary(i, n):
   
    #convertir un int en une liste de 0 et de 1 de longueur au moins n
    #i= a1*2**0+...+an*2**(n-1) s'écrit [a1,...,an] 
    
    binarystr = bin(i) 
    binarylist = []
   
    for j in range(len(binarystr)-1, 1, -1): #
        binarylist.append(int(binarystr[j]))
    
    while len(binarylist) < n :
        binarylist.append(0)
    
    return np.array(binarylist)

def F2_n(n):
    
    #sous-espace vectoriel de R^n dont les vecteurs ne contiennent que des O et des 1
    
    return [tobinary(i,n) for i in range(0,2**n)]

def symetrique(n):
    
    #matrices symétriques de taille nxn contenant que des 0 et des 1
    
    #On commence par obtenir une base
    
    BaseS = []
    
    for i in range(n):
        
        #On ajoute la matrice Sii avec un 1 en position (i,i)
        Sii = np.zeros((n,n)) 
        Sii[(i,i)] = 1 
        BaseS.append(Sii)
        
        for j in range(i+1,n):
            
            # On ajoute la matrice Sij avec un 1 en position (i,j) et un 1 en position (j,i)
            Sij = np.zeros((n,n))
            Sij[(i,j)] = 1
            Sij = Sij + Sij.T
            BaseS.append(Sij)
    
    #Pour S matrice symétrique de taille nxn contenant que des 0 et des 1
    #S = sum[1<=i<=n,i<=j](aij*Sij) où aij dans {0,1}^(len(BaseS))
    Symetriques = [] 
    m = len(BaseS) 
    F2_m = F2_n(m)
    
    for comblin in F2_m:
    
        #On balaie tous les [aij] possibles
    
        S = np.zeros((n,n))
        
        for i in range(m):
        
            S += comblin[i]*BaseS[i] #S = sum[1<=i<=n,i<=j](aij*Sij) 
        
        Symetriques.append(S)
    
    return Symetriques

def sev(A,B):
    
    #sous espace vectoriel engendré par les lignes de la matrice C=(A,B)
    #S = sum[1<=i<=n](ei*ci) où ci est la i-ème ligne de C et (ei) est dans {0,1}^n 
    #on a un ensemble de lignes en sortie
    
    set = []
    taille = np.shape(A)
    n = taille[0]
    F2n = F2_n(n)
    
    C = np.concatenate((A,B),axis=1)
    
    for comblin in F2n:

        #on parcourt les vecteurs (ei)
        
        row = np.zeros(2*n)
        
        for i in range(n):
            
            #on ajoute ei*ci 
            row += comblin[i]*C[i]
        
        set.append(row)
    
    return set
    
    
def pauli_ordre_2(i,j):
    
    #les matrices de Pauli sont codées de manière binaire
    #00 : I
    #10 : X
    #01 : Y
    #11 : Z
    
    if i%2==0 and j%2==0:
        return 'I'
    elif i%2==1 and j%2==0 :
        return 'X'
    elif i%2==0 and j%2==1 :
        return 'Z'
    else :
        return 'Y'
        
def sub_pauli(A,B):
    
    #on traduit le sous-espace vectoriel engendré par les lignes de la matrice (A,B) en matrice de Pauli selon la fonction précédente
    set = sev(A,B)
    taille = np.shape(A)
    n = taille[0]
    paulisubset = []
    
    for row in set:
        
        #pour chaque ligne du sous-espace vectoriel
        #on définit une matrice de pauli à n qubit
        
        pauli_ordre_n = ''
        
        for i in range(n):
            
            #la ième matrice d'ordre 2 du n-qubit est donnée par les valeurs en position i et n+i de la ligne du sous-ensemble
            pauli_ordre_n += pauli_ordre_2(row[i],row[n+i])
        
        paulisubset.append(pauli_ordre_n)
    
    return paulisubset
    
def commutative_pauli_subset(n):
    
    CommutativePauliSubsets = []
   
    Symetrique = symetrique(n)
    
    CommutativePauliSubsets.append(sub_pauli(np.zeros((n,n)),np.eye(n)))
    
    for A in Symetrique:
        CommutativePauliSubsets.append(sub_pauli(np.eye(n),A))
    
    return CommutativePauliSubsets
 
def all_commutative_pauli_subset(n):
    
    CommutativePauliSubsets = []
   
    Symetrique = symetrique(n)
    
    for A in Symetrique:
        for B in Symetrique :
            CommutativePauliSubsets.append(sub_pauli(A,B))
    
    return CommutativePauliSubsets
    
 
def all_matrix_set(n):
    
    Symetrique = symetrique(n)[1:]
    m = len(Symetrique)
    d = 2**n
    
    set = []

    for i in range(m) :
        
        
        subset = [Symetrique[i]]
        
        for j in range(m) :
            
            M = Symetrique[j]
            condition = True
        
            for A in subset :
            
                if np.linalg.det(A-M) == 0 :
               
                    condition = False
        
            if condition :
                subset.append(M)
        
        if len(subset)==d-1:
            subset.append(np.zeros((n,n)))
            set.append(subset)
    
    return set

def MUBs(n):
    
    MUBs = []
   
    SetSetMatrix = all_matrix_set(n)
    
    for SetMatrix in SetSetMatrix :
       
        MUB = []
        MUB.append(sub_pauli(np.zeros((n,n)),np.eye(n)))
    
        for A in SetMatrix:
        
            MUB.append(sub_pauli(np.eye(n),A))
        
        MUBs.append(MUB)
    
    return MUBs
    
def one_matrix_set(n):
    
    Symetrique = symetrique(n)[1:]
    m = len(Symetrique)
    d = 2**n

    for i in range(m) :
        
        subset = [Symetrique[i]]
        
        for j in range(m) :
            
            M = Symetrique[j]
            condition = True
        
            for A in subset :
            
                if np.linalg.det(A-M) == 0 :
               
                    condition = False
        
            if condition :
                subset.append(M)
        
        if len(subset) == d-1:
            subset.append(np.zeros((n,n)))
            return subset

def MUB(n):
    
    SetMatrix = one_matrix_set(n)
    MUB = []
    MUB.append(sub_pauli(np.zeros((n,n)),np.eye(n)))
    
    for A in SetMatrix:
        
        MUB.append(sub_pauli(np.eye(n),A))
    
    return MUB
    