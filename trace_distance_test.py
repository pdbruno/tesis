from qiskit_aer import AerSimulator
from channels import DepolarizingChannel, AmplitudeDampingChannel, MirroredAmplitudeDampingChannel, NoiselessChannel
from averages import haar_measure_average, pauli_eigenvectors_average, average_distance_of_teleportation
from distances import trace_distance
import numpy as np

shots_per_input = 100000
simulator = AerSimulator()

d = np.vectorize(trace_distance, signature='(3),(3)->()')

    
print(f'pauli_eigenvectors_average: {average_distance_of_teleportation(pauli_eigenvectors_average(), d, NoiselessChannel(), NoiselessChannel(), simulator, shots_per_input)}')
print(f'haar_measure_average: {average_distance_of_teleportation(haar_measure_average(1000), d, NoiselessChannel(), NoiselessChannel(), simulator, shots_per_input)}')

