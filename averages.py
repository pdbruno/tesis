from qiskit.quantum_info import random_unitary, DensityMatrix
import numpy as np

zero = DensityMatrix.from_label("0")

def to_bloch(rho: DensityMatrix):
    [[a, b], [c, d]] = rho.data
    return np.array([complex(c + b).real, complex(c - b).imag, complex(a-d).real])

def haar_measure_average(shots):
    for i in range(shots):
        haar_random_operator = random_unitary(2)
        f = lambda qc, psi: qc.append(haar_random_operator.to_instruction(), [psi])
        yield f, to_bloch(zero.evolve(haar_random_operator))
        
def eigenvector_x_mas(qc, psi_alice):
    qc.h(psi_alice)
def eigenvector_x_menos(qc, psi_alice):
    qc.x(psi_alice)
    qc.h(psi_alice)
def eigenvector_y_mas(qc, psi_alice):
    qc.h(psi_alice)
    qc.s(psi_alice)
def eigenvector_y_menos(qc, psi_alice):
    qc.h(psi_alice)
    qc.sdg(psi_alice)
def eigenvector_z_mas(qc, psi_alice):
    pass
def eigenvector_z_menos(qc, psi_alice):
    qc.x(psi_alice)
    
def pauli_eigenvectors_average():
    yield eigenvector_x_mas, np.array([1., 0., 0.])
    yield eigenvector_x_menos, np.array([-1., 0., 0.])
    yield eigenvector_y_mas, np.array([0., 1., 0.])    
    yield eigenvector_y_menos, np.array([0., -1., 0.])    
    yield eigenvector_z_mas, np.array([0., 0., 1.])    
    yield eigenvector_z_menos, np.array([0., 0., -1.])      