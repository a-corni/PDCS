import numpy as np
import sys
import Classic_algebra
import Classic_generators
np.set_printoptions(threshold=sys.maxsize)

def linear_subspace(n):
    
    # Retourne tous les sous-ensembles commutatifs d'opérateurs de Pauli de taille N
    #O(N**2*2**(N**2)+ N**4 + (C(N,2N)*2**N)**2) = O(N**2*2**(N**2) + 2**(6N)) = O(N**2*2**(N**2)) 
    
    # liste des sous-ensembles commutatifs d'opérateurs de Pauli
    linearsubspace = []
    
    # Chaque générateur d'un ensemble commutatif d'opérateur de Pauli s'écrit [A,B] avec A,B des matrices carrées ne contenant que des 0 et des 1
    MnF2 = Mn_F2(n) #O(N**2*2**(N**2))
    
    for A in MnF2 : #n**2 call
       
        for B in MnF2 :#n**2 call
            
            # générateur
            M = np.concatenate((A,B),axis=1) #1 operation 
            
            # on génère un groupe commutatif de matrices de Pauli si M = [A,B] est de rang n
            if np.linalg.matrix_rank(M) == n: # C(2n,n)*2**n matrices
                
                # on regarde si M ne génère pas le même ensemble qu'un générateur déjà observé
                hasappeared = False
                m = len(linearsubspace)
                
                # on parcourt la liste des ensembles déjà générés
                for indexsubspace in range(m): # unknown but maximal C(2n,n)*2**n
                    
                    subspace = linearsubspace[indexsubspace]
                    representant = subspace[0] 
                    
                    # on compare notre générateur à un générateur déjà trouvé
                    if are_equivalent(M, representant): #O(n*2**N)
                    
                        # s'ils génèrent le même groupe
                        # on stocke la matrice M au format (A,B) dans la même liste que le générateur comparé 
                        linearsubspace[indexsubspace].append((A,B))
                        hasappeared = True
                
                if not hasappeared :
                    
                    # sinon
                    # on crée un nouveau groupe de matrices équivalentes. 
                    linearsubspace.append([M, (A,B)])
    
    return linearsubspace
   
def linear_subset_new_condition(N):
    
    #liste des matrices commutatives
    linear_subsets = []
    
    MnF2 = Mn_F2(n) #O(N**2*2**(N**2))
    
    I = np.eye(N,dtype=int)
    zero = np.zeros((N,N),dtype=int)
    
    for A in MnF2 : 
    
        M = np.concatenate((I,A), axis = 1)
        generated_subsets = permutation(M)
            
        for generated_subset in generated_subsets :
            
            # on regarde si M ne génère pas le même ensemble qu'un générateur déjà observé
            hasappeared = False
            m = len(linear_subsets)
            
            # on parcourt la liste des ensembles déjà générés
            for indexsubset in range(m): # unknown but maximal C(2n,n)*2**n
                
                subset = linearsubset[indexsubset]
                representant = subset[0] 
                
                # on compare notre générateur à un générateur déjà trouvé
                if are_equivalent(generated_subset, representant): #O(N*2**N)
                
                    # s'ils génèrent le même groupe
                    # on stocke la matrice M au format (A,B) dans la même liste que le générateur comparé 
                    linearsubspace[indexsubset].append(generated_subset)
                    hasappeared = True
            
            if not hasappeared :
                
                # sinon
                # on crée un nouveau groupe de matrices équivalentes. 
                linearsubspace.append([M])
             
    
    return linear_subsets
    
def linear_subspace_toprint(n):

    fichier = open("linear_subspace.txt", "w")
    linearsubspace = linear_subspace(n)
    understoodable = []
    
    for subspace in linearsubspace:
        
        n = len(subspace)
        (A,B) = subspace[1]
        subspacename = sub_pauli(A,B)
        understoodable.append(subspacename)
        fichier.write(str(subspacename) + "Il y a " + str(n-1)+ " matrices équivalentes :" +"\n")
        
        for i in range(1,n):
            
            understoodable.append(subspace[i])
            fichier.write(str(subspace[i]) + ", "+"\n")
        
        understoodable.append(n-1)
        fichier.write("\n")
    
    fichier.close()
    
    return understoodable