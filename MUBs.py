import numpy as np
import sys
import Classic_algebra
import Classic_generators
import csv
np.set_printoptions(threshold=sys.maxsize)

def all_set_A(N):
    
    # En entrée : S = ensemble des matrices symétriques de taille n
    # On veut construire tous les ensembles de matrices {A_1,...,A_(2^N)} tels que pour 1<=i<j<=N, det(A_i-A_j)!=0
    
    All_sets_A = [] # liste des ensembles de matrices {A_1,...,A_(2^N)}
    
    def set_A_rec(built_A, next_A):
        
        new_S = []
        A_i = built_A[-1]
        for A_j in next_A :
            if np.linalg.det(A_i-A_j)!=0:
                new_S.append(A_j)
        
        for A_j in new_S :
            new_next_A = list(new_S)
            new_built_A = built_A+[A_j]
            
            if len(new_built_A) == 4:
                All_sets_A.append(new_built_A)
            else :
                new_next_A = remove_array_from_list_array(new_next_A, A_j)
                set_A_rec(new_built_A,new_next_A)
            
    S = Symetrique(N)
    
    for A_1 in S:
        
        next_A = list(S)
        next_A = remove_array_from_list_array(next_A, A_1)
        set_A_rec([A_1],next_A)
        
    return All_sets_A
        
def multiple_MUB(N):
    
    MUBs = []
    all_sets_A = all_set_A(N)
    I = np.eye(N,dtype=int)
    zero = np.zeros((N,N),dtype=int)
    
    for set_A in all_sets_A : 
        
        for i in range(0,2*N):
            
            
            MUB_i = []
        
            MUB_i.append(translation_operator(np.concatenate((zero,I),axis=1), i))
        
            for A in set_A :
        
                MUB_i.append(translation_operator(np.concatenate((I,A),axis=1), i)) 
    
            MUBs.append(MUB_i)
    
    return MUBs
   
def multiple_MUB_understoodable(N):
    
    fichier = open("MUB_Pauli_order_"+str(n)+".csv", "wt")
    MUBCSV = csv.writer(fichier,delimiter=";")
    ecrivainCSV.writerow(["MUB index"] + [str(i+1) for i in range(2**N+1)]
    mubs= multiple_MUB(N)
    
    for MUB in mubs:
        
        csvline = []
        
        for commutative_subgroups in MUB:
            
            X = commutative_subgroups[:,:N]
            Z = commutative_subgroups[:,N:]
             
            csvline.append(sub_pauli(X,Z))
        MUBCSV.writerow(csvline)
    
    fichier.close()

