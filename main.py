from concurrent.futures import ThreadPoolExecutor
from qiskit_aer.primitives import SamplerV2
from noisy_quantum_teleportation_benchmarking.channels import (
    DepolarizingChannel,
    AmplitudeDampingChannel,
    MirroredAmplitudeDampingChannel,
    PhaseDampingChannel,
)
from noisy_quantum_teleportation_benchmarking.distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from noisy_quantum_teleportation_benchmarking.experiment import Experiment
from qiskit.compiler import transpile

from noisy_quantum_teleportation_benchmarking.sampler import HaarMeasureSampler, PauliSampler

channels = [
    AmplitudeDampingChannel(),
    MirroredAmplitudeDampingChannel(),
    DepolarizingChannel(),
    PhaseDampingChannel()
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

"""     
    }, """

sampler = SamplerV2(
    default_shots=2000,
    options={
        "backend_options": {
            "executor": ThreadPoolExecutor(max_workers=10),
            "max_job_size": 1,
            "max_parallel_experiments": 0,
        }}
)

#input_sampler = PauliSampler()
input_sampler = HaarMeasureSampler(500)


transpiler = lambda qc: transpile(qc, sampler._backend)
experiment = Experiment(channel_combinations, exploration_space, input_sampler, transpiler)
experiment.run_with_sampler(sampler, distances, sampler._default_shots, 'haar-avg-sim-2.csv')