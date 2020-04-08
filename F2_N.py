import numpy as np
import sys
import Classic_algebra
import Classic_generators
np.set_printoptions(threshold=sys.maxsize)

def linear_subspace(n):
    
    #O(N**2*2**(N**2)) complexity
    
    linearsubspace = []
    MnF2 = Mn_F2(n) #O(N**2*2**(N**2))
    
    for A in MnF2 : #n**2 call
       
        for B in MnF2 :#n**2 call
    
            M = np.concatenate((A,B),axis=1) #1 operation 
    
            if np.linalg.matrix_rank(M) == n: # C(2n,n)*2**n matrices
                
                hasappeared = False
                m = len(linearsubspace) #1 operation
                
                for indexsubspace in range(m): #unknown but maximal C(2n,n)*2**n
                    
                    subspace = linearsubspace[indexsubspace]
                    representant = subspace[0]
                    
                    if are_equivalent(M, representant): #O(n**2*2**N)
                        linearsubspace[indexsubspace].append((A,B))
                        hasappeared = True
                
                if not hasappeared :
                    linearsubspace.append([M, (A,B)])
    
    return linearsubspace
   
def linear_subspace_understoodable(n):
    
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