import numpy as np
from math import log
import sys
np.set_printoptions(threshold=sys.maxsize)


def tobinary(i, n):
    
   
    # convertir un int en une liste de '0' et de '1' de longueur n (attention se sont des str) 
    # On s'imposera i < 2**n
    # i= a1*2**0+...+an*2**(n-1) s'écrit [an,...,a1] 
    #O(1)
    
    #fonction qui convertit un entier en une chaîne de caractère représentant sa décomposition en base 2
    length = '{0:0'+str(n)+'b}'
    binary_string = length.format(i)
    return list(binary_string)

    
def F2_n(n):
    
    #sous-espace vectoriel de R^n dont les vecteurs ne contiennent que des O et des 1
    #O(2**n)
    return [tobinary(i,n) for i in range(0,2**n)]

def Mn_F2(n):
    
    #Matrices de taille NxN composées seulement de 0 et de 1
    #Complexité en O(2**(N**2)*N**2)
    
    #On construit une base de matrices de Mn
    #Se sont les matrices dont les composantes sont nulles partout sauf en une position (i,j).
    Base = []
    
    for i in range(n): # n calls
        
        for j in range(n): # n calls
            
            
            M = np.zeros((n,n),int)
            M[i,j]=1 #O(1)
            Base.append(M) #O(1)
    
    #A partir de cette base on construit les matrices de taille NxN composées seulement de 0 et de 1 
    #On fait des combinaisons linéraires avec que des 0 et des 1 des matrices de la base       
    MnF2 = []
    F2n2 = F2_n(n**2) #listes des combinaisons linéaires O(2**(N**2))
    
    for comblin in F2n2: #2**(N**2) calls
        
        M = np.zeros((n,n),int)
        
        for i in range(n**2): #N**2 calls
            
            M += int(comblin[i])*Base[i] #O(1)
    
        MnF2.append(M) #O(1)
    
    return MnF2

def Symetrique(n):
    
    #matrices symétriques de taille nxn contenant que des 0 et des 1
    #On commence par obtenir une base
    #O(n(n+1)/2*2**(n(n+1)/2))
    
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
    F2_m = F2_n(m) #O(2**(n(n+1)/2))
    
    for comblin in F2_m: #2**(n(n+1)/2) calls
    
        #On balaie tous les [aij] possibles
    
        S = np.zeros((n,n),int)
        
        for i in range(m): #n(n+1)/2
        
            S += int(comblin[i])*BaseS[i] #S = sum[1<=i<=n,i<=j](aij*Sij) 
        
        Symetriques.append(S)
    
    return Symetriques

def translation_operator(M, i):
    
    # retourne la matrice N telle que N[k,j] = M[k, j+i[n]]
    
    taille = np.shape(M)
    n = 2*taille[0]
    
    return np.concatenate((M[:,n-i:],M[:,:n-i]), axis = 1)

def remove_array_from_list_array(list, test_array):
    
    return [ M for M in list if not (M==test_array).all()]