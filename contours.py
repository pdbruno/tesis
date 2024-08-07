from qiskit_aer import AerSimulator
from channels import DepolarizingChannel, AmplitudeDampingChannel, MirroredAmplitudeDampingChannel
from averages import pauli_eigenvectors_average, average_distance_of_teleportation
from distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from concurrent.futures import ThreadPoolExecutor

shots_per_input = 200000
simulator = AerSimulator()
exc = ThreadPoolExecutor(max_workers=10)
simulator.set_options(executor=exc)
simulator.set_options(max_job_size=1)


channels = [AmplitudeDampingChannel, MirroredAmplitudeDampingChannel, DepolarizingChannel]
distances = list(map(lambda d: np.vectorize(d, signature='(3),(3)->()'), [trace_distance, fidelity, affinity, wooters_distance]))
    
ps = np.linspace(0, 1, 11)

avg = pauli_eigenvectors_average()

with open('res.txt', '+w') as f:
    for chA in channels:
        for chB in channels:
            for d in distances:
                f.write('--------------')
                f.write('\n')
                f.write(chA.label)
                f.write('\n')
                f.write(chB.label)
                f.write('\n')
                f.write(d.__name__)
                f.write('\n')
                values = [(x.item(), y.item(), average_distance_of_teleportation(avg, d, chA(x), chB(y), simulator, shots_per_input).item()) for x in ps for y in ps]
                f.write(str(values))
                f.write('\n')
