import numpy as np
import sys
import Classic_algebra
import Classic_generators
np.set_printoptions(threshold=sys.maxsize)

def conjecture_linear_subspace(n):
    
    #O(N(N+1)/2**(N(N+1)/2)) complexity
    
    linearsubspace = []
    symetrique = Symetrique(n) #O(N(N+1)/2*2**(N(N+1)/2))
    
    for S in symetrique : #n(n+1)/2 calls
        
        I = np.eye(n,dtype=int)
        M = np.concatenate((I,S),axis=1) 
            
        for i in range(1,2*n+1): # 2*n calls
                
                M_i = translation_operator(M,i) #O(1)

                hasappeared = False
                m = len(linearsubspace) 
                
                for indexsubspace in range(m): #unknown but maximal 2*n*2**(n(n+1)/2)
                    
                    subspace = linearsubspace[indexsubspace]
                    representant = subspace[0]
                    
                    if are_equivalent(M_i, representant): #O(n**2*2**N)
                        linearsubspace[indexsubspace].append(M_i)
                        hasappeared = True
                
                if not hasappeared :
                    linearsubspace.append([M_i, M_i])
    
    return linearsubspace
   
def conjecture_linear_subspace_understoodable(n):
    
    fichier = open("conjecture_linear_subspace"+str(n)+".txt", "w")
    linearsubspace = conjecture_linear_subspace(n)
    understoodable = []
    
    fichier.write("Il y a "+ str(len(linearsubspace))+" groupes commutatifs de matrices de Pauli à " +str(n) +"-qubits"+ "\n")
    
    for subspace in linearsubspace:
        
        m = len(subspace)
        M = subspace[1]
        M1 = M[:,:n]
        M2 = M[:,n:]
        subspacename = sub_pauli(M1,M2)
        understoodable.append(subspacename)
        fichier.write(str(subspacename) + "Il y a " + str(m-1)+ " matrices équivalentes :" +"\n")
        
        for i in range(1,m):
            
            understoodable.append(subspace[i])
            fichier.write(str(subspace[i]) + ", "+"\n")
        
        understoodable.append(m-1)
        fichier.write("\n")
    
    fichier.close()
    
    return understoodable

