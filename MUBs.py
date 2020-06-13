import numpy as np
import time
import sys
import math
from Classic_algebra import *
from Classic_generators import *
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
            
            if len(new_built_A) == 2**N:
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
    start = time.time()
    
    for set_A in all_sets_A : 
        
        # On a {A_1,...,A_(2^N)} tels que pour 1<=i<j<=N, det(A_i-A_j)!=0
        
        computing_time = time.time()-start
        
        if computing_time < 180 :
            
            for i in range(0,2*N):
                
                # MUB_i = {T^i([1|A_1]),...,T^i([1|A_(2^N)])} 
                MUB_i = []
            
                MUB_i.append(translation_operator(np.concatenate((zero,I),axis=1), i))
            
                for A in set_A :
            
                    MUB_i.append(translation_operator(np.concatenate((I,A),axis=1), i)) 
                
                #on ne la rajoute à la liste que si on n'a pas déjà réalisé cette répartition
                
                already_counted = False
                
                for MUB in MUBs :
                    
                    same_repartition = True
                    
                    for commutative_subgrp_MUB_i in MUB_i:
                        
                        is_in_MUB = False
                        
                        for commutative_subgrp_MUB in MUB :
                        
                            if are_equivalent(commutative_subgrp_MUB_i, commutative_subgrp_MUB):
                                is_in_MUB = True
                                break
                        
                        if is_in_MUB == False :
                            same_repartition = False
                            break
                    
                    if same_repartition :
                        already_counted = True
                        break
                    
                if not already_counted :
                        
                    MUBs.append(MUB_i)
            
    return MUBs
   
def multiple_MUB_toprint(N):
    
    fichier = open("MUB_Pauli_order_"+str(N)+".csv", "wt")
    MUBCSV = csv.writer(fichier,delimiter=";")
    MUBCSV.writerow(["index "+ str(i+1) for i in range(2**N+1)])
    MUBs= multiple_MUB(N)
    
    for MUB in MUBs:
        
        csvline = []
        
        for commutative_subgroups in MUB:
            
            X = commutative_subgroups[:,:N]
            Z = commutative_subgroups[:,N:]
             
            csvline.append(sub_pauli(X,Z))
        MUBCSV.writerow(csvline)
    
    fichier.close()

