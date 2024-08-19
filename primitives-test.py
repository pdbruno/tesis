from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2
from channels import (
    DepolarizingChannel,
    AmplitudeDampingChannel,
    MirroredAmplitudeDampingChannel,
    NoiselessChannel,
)
from averages import get_score_protocol_distance, pauli_eigenvectors_average, to_bloch
from distances import affinity, trace_distance, wooters_distance, fidelity
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from pub_generator import generate_pub
from qiskit.compiler import transpile
import pandas as pd
from quantum_state_tomography import get_probabilities_and_states
from typing import cast


def get_counts(array):
    bitstrings, counts = np.unique(
        array,
        axis=0,
        return_counts=True,
    )
    return {
        tuple(bitstring.data): count.item()
        for bitstring, count in zip(bitstrings, counts)
    }


""" channels = [
    AmplitudeDampingChannel(),
    MirroredAmplitudeDampingChannel(),
    DepolarizingChannel(),
] """

channels = [AmplitudeDampingChannel()]
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
    default_shots=10000,
)
transpiler = lambda qc: transpile(qc, sampler._backend)
avg = pauli_eigenvectors_average()
pubs = [
    generate_pub(chA, chB, exploration_space, avg, transpiler)
    for chA, chB in channel_combinations
]

results = sampler.run(pubs).result()


def average_distance_of_teleportation(formatted_results, rho_A):

    p_is, rho_B_is = get_probabilities_and_states(
        [get_counts(formatted_results[idx_basis].T) for idx_basis in range(3)],
        sampler._default_shots,
    )
    return [get_score_protocol_distance(p_is, rho_A, rho_B_is, d) for d in distances]


for pub_result, (chA, chB) in zip(results, channel_combinations):
    length = pub_result.data["input_meas"].shape[0]
    formatted_results = np.split(
        np.array(
            np.split(
                np.swapaxes(
                    np.array(
                        [
                            np.squeeze(pub_result.data["input_meas"].array),
                            np.squeeze(pub_result.data["alice_meas"].array),
                            np.squeeze(pub_result.data["bob_meas"].array),
                        ]
                    ),
                    0,
                    1,
                ),
                cast(int, length / 3),
            )
        ),
        cast(int, length / 3 / 6),
    )
    # (121, 6, 3, 3, 1000) = (noise_grid, input_sample, measurement_basis, measured_qubits, shots)

    df = pd.DataFrame(
        {
            "chA": [],
            "chB": [],
            "pA": [],
            "pB": [],
            "d": [],
            "score": [],
        }
    )

    for res_for_noise_configuration, (pA, pB) in zip(
        formatted_results, exploration_space
    ):
        scores = np.array(
            [
                average_distance_of_teleportation(
                    res_for_input_state, to_bloch(phi_eulerian_angles)
                )
                for res_for_input_state, phi_eulerian_angles in zip(
                    res_for_noise_configuration, avg
                )
            ]
        ).mean(axis=0)
        df = pd.concat(
            [
                pd.DataFrame(
                    [
                        [chA.label, chB.label, pA, pB, d.__name__, score]
                        for d, score in zip(distances, scores)
                    ],
                    columns=df.columns,
                ),
                df,
            ],
            ignore_index=True,
        )
    df.to_csv('adc-adc-meas.csv')