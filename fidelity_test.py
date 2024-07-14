from qiskit_aer import AerSimulator
from channels import depolarizing_channel, amplitude_damping_channel, mirrored_amplitude_damping_channel
from averages import haar_measure_average, pauli_eigenvectors_average, average_distance_of_teleportation
from distances import fidelity
import numpy as np

def sin_ruido(qc, qubit):
    pass


shots_per_input = 10000
simulator = AerSimulator()

d = np.vectorize(fidelity, signature='(3),(3)->()')

    
print(f'pauli_eigenvectors_average: {average_distance_of_teleportation(pauli_eigenvectors_average(), d, sin_ruido, sin_ruido, simulator, shots_per_input)}')
print(f'haar_measure_average: {average_distance_of_teleportation(haar_measure_average(100), d, sin_ruido, sin_ruido, simulator, shots_per_input)}')

