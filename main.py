from qiskit_aer.primitives import SamplerV2
from channels import (
    DepolarizingChannel,
    AmplitudeDampingChannel,
    MirroredAmplitudeDampingChannel,
)
from averages import pauli_eigenvectors_average
from distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from experiment import Experiment
from qiskit.compiler import transpile

channels = [
    AmplitudeDampingChannel(),
    MirroredAmplitudeDampingChannel(),
    DepolarizingChannel(),
]
channel_combinations = [(chA, chB) for chA in channels for chB in channels]
distances = list(
    map(
        lambda d: np.vectorize(d, signature="(3),(3)->()"),
        [trace_distance, fidelity, affinity, wooters_distance],
    )
)
ps = np.linspace(0, 1, 11)
exploration_space = [(x, y) for x in ps for y in ps]

"""     options={
        "backend_options": {
            "executor": ThreadPoolExecutor(max_workers=10),
            "max_job_size": 1,
            "max_parallel_experiments": 0,
        }
    }, """

sampler = SamplerV2(
    default_shots=2000,
)
input_sampler = pauli_eigenvectors_average()
transpiler = lambda qc: transpile(qc, sampler._backend)
experiment = Experiment(channel_combinations, exploration_space, input_sampler, transpiler)
df = experiment.run_with_sampler(sampler, distances, sampler._default_shots)
df.to_csv('simulation_pauli_eigenvectors_average.csv')