import numpy as np
from cmath import *
import csv
from Classic_algebra import *

# Number of qubits

N = 3

# N-qubit gate U

U = np.zeros([2**N,2**N], complex)
for i in range (2**N):
    for j in range(2**N):
        U[i,j] = e**(2*pi*1j*i*j/2**N)/2**(N/2)

print("Your matrix is :")
print(U)

# Extract H matrix from gate

H = extract_hermitian_from_gate(U)

# All N-qubit Pauli operators

(pauli_names, pauli_matrices) = SU_2_N(N)

# Decomposition of H over the set of N-qubit operators

decomposition = [HS(H,pauli_matrix,N) for pauli_matrix in pauli_matrices]

#Showing results
#Showing [[O_i, alpha_i]]

decomposition_output = order(fusion_list(pauli_names, decomposition))

#Analysis

tokeep = []
n = len(decomposition_output)

for i in range(n-1,-1,-1):
    x = decomposition_output[i]
    if abs(x[1])> 10**(-10):
        y = [x[0], x[1].real]
        tokeep.append(y)
    else :
        break

print(tokeep)

#MUB inspired classification :

fichier = open("MUB_Pauli_order_"+str(N)+".csv", "rt")
MUBCSV = csv.reader(fichier,delimiter=";")
possible_partitioning = { i:[] for i in range(2**N+2)}
i = 0
for row in MUBCSV:
    if i!=0:
        partitioning = []
        for base in row :
            commutative_set = []
            for x in tokeep : 
                if x[0] in base:
                    commutative_set.append(x)
            partitioning.append(commutative_set)
        k = len(partitioning)
        possible_partitioning[k].append(partitioning)
    else :
        i+=1

best_index = 0
for k in range(1,2**N+2):
    if len(possible_partitioning[k])!=0:
        best_index = k
        break

print("Les solutions permettant le partitionnement à plus faible coût sont :"+"\n")

index_solution = 0

for solution in possible_partitioning[best_index]:
    
    index_solution+=1
    print("Solution "+str(index_solution)+" : "+"\n")
    index_gate = 0
    print("Phase factor :")
    print(solution[0][0])
    for gate in solution :
        if len(gate)!=1:
            index_gate+=1
            print("Décomposition Hamiltonien porte "+str(index_gate) + " : ")
            print(gate[1:])

#First Iterative Classification :

cover = [] 
phase_factor = []
identity = ''
for _ in range(N):
    identity+='I'

for pauli_operator_phase in tokeep :
    
    (pauli_operator,phase) = (pauli_operator_phase[0],pauli_operator_phase[1])
    
    if pauli_operator == identity:
        
        phase_factor.append([pauli_operator_phase])
    
    else :
        n = len(cover)
        fitting_subsets = []
        for i in range(n) :
            subset = cover[i]
            commute_with_all = True
            for operator in subset:
                
                if not commute(pauli_operator, operator[0], N) :
                    commute_with_all = False
            if commute_with_all :
                fitting_subsets.append(i)
            
        m = len(fitting_subsets)
        if m == 0 :
            cover.append([pauli_operator_phase])
        else :
            for index in fitting_subsets:
                cover[index].append([pauli_operator, phase/m])
print("Decomposition using an iterative cover")
print(phase_factor+cover)
print("Number of gates :")
print(len(phase_factor+cover))

#Second Iterative Classification :

cover = [] 
phase_factor = []
identity = ''
for _ in range(N):
    identity+='I'

for pauli_operator_phase in tokeep :
    
    (pauli_operator,phase) = (pauli_operator_phase[0],pauli_operator_phase[1])
    
    if pauli_operator == identity:
        
        phase_factor.append([pauli_operator_phase])
    
    else :
        n = len(cover)
        fitting_subsets = []
        for i in range(n) :
            subset = cover[i]
            commute_with_all = True
            for operator in subset:
                
                if not commute(pauli_operator, operator[0], N) :
                    commute_with_all = False
            if commute_with_all :
                fitting_subsets.append(i)
            
        m = len(fitting_subsets)
        cover.append([pauli_operator_phase])
        for index in fitting_subsets:
                cover[index].append([pauli_operator, phase])

for pauli_operator_phase in tokeep :
    
    (pauli_operator,phase) = (pauli_operator_phase[0],pauli_operator_phase[1])
    
    if not pauli_operator == identity:
    
        n = len(cover)
        
        for i in range(n) :
            
            subset = cover[i]
            commute_with_all = True
            
            if pauli_operator_phase in subset :
            
                break
            
            else :
                    
                for operator in subset:
                    
                    if not commute(pauli_operator, operator[0], N) :
               
                        commute_with_all = False
               
                if commute_with_all :
               
                    cover[i].append([pauli_operator, phase])            

Final_cover = []

while len(cover)!=0:
    
    weights = [norm(subset) for subset in cover]
    m = max(weights)
    n = len(weights)
    for i in range(n):
        if weights[i] == m:
            best = i
            break
    best_subset = cover[i]
    Final_cover.append(list(best_subset))
    n = len(cover)
    for i in range(n-1, -1, -1):
        subset = cover[i]
        m = len(subset)
        todel= []
        for j in range(m):
            if subset[j] in best_subset :
                todel.insert(0, j)
                
        for j in todel:
            del subset[j]
            
        if len(subset)==0:
            del cover[i]
print("Final_cover :")
print(phase_factor+Final_cover)
print("Number of gates :")
print(len(phase_factor+Final_cover))  
     

    
    
