from channels import (
    DepolarizingChannel,
    AmplitudeDampingChannel,
    MirroredAmplitudeDampingChannel,
    PhaseDampingChannel,
)
from distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from experiment import Experiment
from qiskit.compiler import transpile
from sampler import HaarMeasureSampler, PauliSampler
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

service = QiskitRuntimeService()
backend = service.least_busy(simulator=False, operational=True)

sampler = Sampler(mode=backend)

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

input_sampler = PauliSampler()
#input_sampler = HaarMeasureSampler(500)


transpiler = lambda qc: transpile(qc, backend)
experiment = Experiment(channel_combinations, exploration_space, input_sampler, transpiler)
experiment.run_with_sampler(sampler, distances, 2000, 'ADC-pauli-quantum.csv')