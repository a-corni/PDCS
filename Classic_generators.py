import Classic_algebra
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

def sev(A,B):
    
    #sous espace vectoriel engendré par les lignes de la matrice C=(A,B)
    #S = sum[1<=i<=n](ei*ci) où ci est la i-ème ligne de C et (ei) est dans {0,1}^n 
    #on a un ensemble de lignes en sortie
    #O(n*2**n) complexity
    
    set = []
    taille = np.shape(A)
    n = taille[0]
    F2n = F2_n(n) #O(2**n)
    
    C = np.concatenate((A,B),axis=1)
    
    for comblin in F2n: 
        
        #on parcourt les vecteurs (ei)
        #2**n calls
        
        row = np.zeros(2*n, int)
        
        for i in range(n): #n calls
            
            #on ajoute ei*ci
            #n calls 
            row += int(comblin[i])*C[i]
        
        set.append(row)
    
    set = np.remainder(set,2)
    
    return np.array(set)



def are_equivalent(M, N):
    
    # M, N are two matrices of Mnx2n(F2)
    # We want to know if they generate the same subspace
    # This happens if all the lines of N belong to the subspace generated by M and of the lines of M belong to the subspace generated by N
    # O(N*2**N) complexity
    
    taille = np.shape(M)
    n = taille[0]
    M1 = M[:,:n]
    M2 = M[:,n:]
    SM = sev(M1,M2) #O(n*2**n)
    
    
    for generatorsN in N: #N calls
       
        is_elementM = False
       
        for elements in SM : #2**N calls
        
            if np.array_equal(generatorsN,elements): #O(1)
               
                is_elementM = True 
       
        if not is_elementM : 
            
            return False
           
    N1 = N[:,:n]
    N2 = N[:,n:]
    SN = sev(N1,N2) #O(N*2**N)

    for generatorsM in M: #N call
        
        is_elementN = False
        
        for elements in SN : #2**N calls
        
            if np.array_equal(generatorsM,elements): #O(1)
        
                is_elementN=True
        
        if not is_elementN :
        
            return False
    
    return True

def pauli_ordre_2(i,j):
    
    #les matrices de Pauli sont codées de manière binaire
    #00 : I
    #10 : X
    #01 : Y
    #11 : Z
    #O(1) complexity
    
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
    #O(N*2**N) complexity
    
    set = sev(A,B) # O(N*2**N)
    taille = np.shape(A)
    n = taille[0]
    paulisubset = []
    
    for row in set: #N calls
        
        #pour chaque ligne du sous-espace vectoriel
        #on définit une matrice de pauli à n qubit
        
        pauli_ordre_n = ''
        
        for i in range(n): #N calls
            
            #la ième matrice d'ordre 2 du n-qubit est donnée par les valeurs en position i et n+i de la ligne du sous-ensemble
            pauli_ordre_n += pauli_ordre_2(row[i],row[n+i])
        
        paulisubset.append(pauli_ordre_n)
    
    return paulisubset
 