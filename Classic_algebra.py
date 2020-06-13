import numpy as np
from math import log
import sys
from itertools import permutations
from cmath import *
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

def permutation(M):
    
   N = np.transpose(M)
   all_permutation = []
   permute_columns = permutations(N)
   
   for tuple in permute_columns :
       
        matrix = []
       
        for column in tuple :

            matrix.append(column)
        
        final_matrix_transposed = np.array(matrix)
        all_permutation.append(np.transpose(final_matrix_transposed))
        
   return all_permutation
   
def tensorprod(matrices):
    
    finalmatrix = matrices.pop(0)
    
    for matrix in matrices :
        
        finalmatrix = np.kron(finalmatrix, matrix)
    
    return finalmatrix
    
def generatematrix(eigenvalues, eigenvectors):
    
    P = np.array(eigenvectors)
    D = np.diag(eigenvalues)
    invP = np.linalg.inv(P)
    
    left = np.matmul(P,D)
    
    return np.matmul(left, invP)
    
def generategate(eigenvaluesH, eigenvectors):
    
    #mwesh ca marche pas bien
    
    P = np.array(eigenvectors)
    eigenvaluesU = [exp(x*1j) for x in eigenvaluesH]
    D = np.diag(eigenvaluesU)
    invP = np.linalg.inv(P)
    
    left = np.matmul(P,D)
    
    return np.matmul(left, invP)
    
def extract_hermitian_from_gate(U):
    
    (eigenvaluesU, P) = np.linalg.eig(U)
    eigenvaluesH = [phase(x) for x in eigenvaluesU]
    
    return generatematrix(eigenvaluesH,P)

def HS(a,b,N):
    
    return np.trace(np.matmul(a,np.conj(b).T))/2**N

def SU_2():
    return (['I','X','Y','Z'], [np.array([[1,0],[0,1]]), np.array([[0, 1],[1, 0]]), np.array([[0, -1*1j],[1j, 0]]), np.array([[1,0],[0,-1]])])

def SU_2_N(N):
    
    (SU_2_names,SU_2_matrices) = SU_2()
    Paulimatrices = list(SU_2_matrices)
    Paulinames = list(SU_2_names)
    
    for i in range(N-1):
        
        m = len(Paulimatrices)
        
        for j in range(m):
            
            Paulimatrix = Paulimatrices.pop(0)
            Pauliname = Paulinames.pop(0)
            
            for k in range(4):
            
                Paulimatrices.append(tensorprod([Paulimatrix, SU_2_matrices[k]]))
                Paulinames.append(Pauliname+SU_2_names[k])
    
    return (Paulinames, Paulimatrices)
    

def fusion_list(L1, L2):
    
    M = []
    n = len(L1)
    
    for i in range(n):
        
        M.append([L1[i],L2[i]])
        
    return M
    
def order(list):
    
    ordered_list = [list[0]]
    n = len(list)
    
    for i in range(1,n):
        
        j = 0
        m = len(ordered_list)
        x = abs(list[i][1])
        
        while j<m and x > abs(ordered_list[j][1]):
                
            j+=1    
            
        ordered_list.insert(j,list[i])
   
    return ordered_list

def eig_projectors(H,N):
    
    (eigenvaluesH, psi) = np.linalg.eig(H)
    
    E = []
    
    for i in range(2**N):
        
        E_i = np.zeros((2**N,2**N))
        psi_i = psi[i]
        
        for j in range(2**N):
            
            E_i[:,j] = psi_i[j]*psi_i
        
        E.append(E_i)
    
    return (eigenvaluesH, E)

def epsilon(N):
    
    eps = np.zeros((N,2**N),int)
    for k in range (1,N+1):
        for i in range(2**N):
            if (i)%(2**(N-k+1))<2**(N-k):
                eps[k-1,i] = 1
            else :
                eps[k-1,i] = -1
    return eps

def alpha(N):
    
    alpha = np.zeros((2**N,2**N),int)
    eps = epsilon(N)
    F2 = F2_n(N)
    
    for k in range(2**N):
        for i in range(2**N):
            aki = 1
            for j in range(N):
                aki*=eps[j,i]**(int(F2[k][j]))
            alpha[k,i] = aki
   
    return alpha
    
def prod_unit(A, B):
    
    if A == 'I':
        return B
    elif B =='I' :
        return A
    elif A == B:
        return 'I'
    elif ((A in ['X','Y']) and (B in ['X','Y'])):
        return 'Z'
    elif ((A in ['X','Z']) and (B in ['X','Z'])):
        return 'Y'
    elif ((A in ['Y','Z']) and (B in ['Y','Z'])):
        return 'X'

def prod(A,B,N):
    
    product = ''
    for i in range(n):
        product+=prod_unit(A[i],B[i])
    return product
        
def commute(A1, A2, n):
    
    flip = 0
    for i in range(n):
        A = A1[i]
        B = A2[i]
        if not (A == B or A == 'I' or B =='I'):
            flip +=1
    if flip%2==0:
        return True
    return False
    
def norm(subset):
    
    S = 0
    
    for x in subset :
        
        S += x[1]**2
    
    return S
