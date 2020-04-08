import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)


def tobinary(i, n):
   
    #convertir un int en une liste de 0 et de 1 de longueur n
    #i= a1*2**0+...+an*2**(n-1) s'écrit [a1,...,an] 
    #O(n+i) complexity
    
    binarystr = bin(i) #O(i)
    binarylist = []
   
    for j in range(len(binarystr)-1, 1, -1): #i calls
        binarylist.append(int(binarystr[j]))#O(1)
    
    while len(binarylist)!=n : #n-i calls
        binarylist.append(0) #O(1)
    
    return np.array(binarylist)

def F2_n(n):
    
    #sous-espace vectoriel de R^n dont les vecteurs ne contiennent que des O et des 1
    #O(n*2**n)
    return [tobinary(i,n) for i in range(0,2**n)]

def Mn_F2(n):
    
    #Matrices de taille NxN composées seulement de 0 et de 1
    #Complexité en O(2**(N**2)*N**2)
    
    Base = []
    
    for i in range(n): # n calls
        
        for j in range(n): #n calls
        
            M = np.zeros((n,n),int)
            M[i,j]=1 #O(1)
            Base.append(M) #O(1)
            
    MnF2 = []
    F2n2 = F2_n(n**2)
    
    for comblin in F2n2: #2**(N**2) calls
        
        M = np.zeros((n,n),int)
        
        for i in range(n**2): #N**2 calls
            
            M += comblin[i]*Base[i] #O(1)
    
        MnF2.append(M) #O(1)
    
    return MnF2

def Symetrique(n):
    
    #matrices symétriques de taille nxn contenant que des 0 et des 1
    #On commence par obtenir une base
    #O(n(n+1)/2*2**(n(n+1)/2)
    
    BaseS = []
    
    for i in range(n): 
        
        # n calls
        #On ajoute la matrice Sii avec un 1 en position (i,i)
        Sii = np.zeros((n,n),int) 
        Sii[(i,i)] = 1 #O(1)
        BaseS.append(Sii)
        
        for j in range(i+1,n): 
        
            #n-i-1 calls
            # On ajoute la matrice Sij avec un 1 en position (i,j) et un 1 en position (j,i)
            Sij = np.zeros((n,n),int)
            Sij[(i,j)] = 1 #O(1)
            Sij[(j,i)] = 1 #O(1)
            BaseS.append(Sij)
    
    #Pour S matrice symétrique de taille nxn contenant que des 0 et des 1
    #S = sum[1<=i<=n,i<=j](aij*Sij) où aij dans {0,1}^(len(BaseS))
    Symetriques = [] 
    m = len(BaseS) #m = n(n+1)/2 
    F2_m = F2_n(m) #O(n(n+1)/2*2**(n(n+1)/2))
    
    for comblin in F2_m: #n(n+1)/2 calls
    
        #On balaie tous les [aij] possibles
    
        S = np.zeros((n,n),int)
        
        for i in range(m):
        
            S += comblin[i]*BaseS[i] #S = sum[1<=i<=n,i<=j](aij*Sij) 
        
        Symetriques.append(S)
    
    return Symetriques

def translation_operator(M, i):
    
    taille = np.shape(M)
    n = 2*taille[0]
    
    return np.concatenate((M[:,n-i:],M[:,:n-i]), axis = 1)

def remove_array_from_list_array(list, test_array):
    
    return [ M for M in list if not (M==test_array).all()]