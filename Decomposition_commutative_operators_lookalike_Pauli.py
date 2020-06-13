import numpy as np
from cmath import *

from Classic_algebra import *

# Number of qubits

N = 2

# N-qubit gate U

U = np.zeros([2**N,2**N], complex)
for i in range (2**N):
    for j in range(2**N):
        U[i,j] = e**(2*pi*1j*i*j/2**N)/2**(N/2)

print("Your matrix is :")
print(U)

# Extract H matrix from gate

H = extract_hermitian_from_gate(U)
(eig,E) = eig_projectors(H,N)

#Computation of the commutative subset of matrices equivalent to Pauli operators of N-qubit
A = alpha(N)
Commutative_subset = []

for k in range(2**N):
    
    O_k = np.zeros((2**N,2**N))
    
    for i in range(2**N):
        
        O_k += A[k,i]*E[i]
    
    Commutative_subset.append(O_k)

#Computation of the projection of H over each matrix O
decomposition_H = []

for i in range(2**N):
   
    coeff = 0.0
   
    for k in range(2**N):
        coeff += eig[k]*A[k,i]
    decomposition_H.append(coeff/2**N)

print("Decomposition of H :")
print(fusion_list(Commutative_subset, decomposition_H))

print("Decomposition of U :")
print([np.exp(decomposition_H[i]*Commutative_subset[i]) for i in range(1,2**N)])