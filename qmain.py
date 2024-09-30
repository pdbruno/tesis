from noisy_quantum_teleportation_benchmarking.QPU_optimized.experiment_QPU_optimized import ExperimentQPUOptimized
from noisy_quantum_teleportation_benchmarking.channels import (
    DepolarizingChannel,
    AmplitudeDampingChannel,
    MirroredAmplitudeDampingChannel,
    PhaseDampingChannel,
)
from noisy_quantum_teleportation_benchmarking.distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from qiskit.compiler import transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

service = QiskitRuntimeService()
backend = service.least_busy(simulator=False, operational=True)

sampler = Sampler(mode=backend) #type: ignore
# Turn on gate twirling. Requires qiskit_ibm_runtime 0.23.0 or later.
sampler.options.twirling.enable_gates = True #type: ignore

sampler.options.dynamical_decoupling.enable = True #type: ignore
sampler.options.dynamical_decoupling.sequence_type = "XpXm" #type: ignore

adc = AmplitudeDampingChannel()
channel_combinations = [(adc, adc)]
distances = list(
    map(
        lambda d: np.vectorize(d, signature="(3),(3)->()"),
        [fidelity],
    )
)
ps = np.linspace(0, 1, 11)
exploration_space = [(x, x) for x in ps]

transpiler = lambda qc: transpile(qc, backend, optimization_level=3)
experiment = ExperimentQPUOptimized(channel_combinations, exploration_space, transpiler)
experiment.run_with_sampler(sampler, distances, 2000, 'ADC-pauli-quantum-QPU-optimized-EM-optim3.csv')