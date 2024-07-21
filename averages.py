from typing import Sequence
from qiskit.quantum_info import random_unitary, DensityMatrix, Operator
import numpy as np
from channels import BaseChannel
from teleportation_circuit import run_simulation
from quantum_state_tomography import get_probabilities_and_states
from qiskit.providers import BackendV2

zero = DensityMatrix.from_label("0")


def to_bloch(rho: DensityMatrix):
    return np.array(
        [
            (rho.data[1, 0] + rho.data[0, 1]).real,
            (rho.data[1, 0] - rho.data[0, 1]).imag,
            (rho.data[0, 0] - rho.data[1, 1]).real,
        ]
    )


def haar_measure_average(shots) -> Sequence[Operator]:
    return (random_unitary(2) for _ in range(shots)) # type: ignore


eigenvector_x_mas = Operator.from_label("H")
eigenvector_x_menos = Operator.from_label("X") & Operator.from_label("H")
eigenvector_y_mas = Operator.from_label("H") & Operator.from_label("S")
eigenvector_y_menos = Operator.from_label("H") & Operator.from_label("S").adjoint()
eigenvector_z_mas = Operator.from_label("I")
eigenvector_z_menos = Operator.from_label("X")


def pauli_eigenvectors_average():
    return [
        eigenvector_x_mas,
        eigenvector_x_menos,
        eigenvector_y_mas,
        eigenvector_y_menos,
        eigenvector_z_mas,
        eigenvector_z_menos,
    ]


def get_score_protocol_distance(p_is, rho_a, rho_B_is, distance):
    return p_is @ distance(rho_a, rho_B_is)


def average_distance_of_teleportation(
    input_sampler: Sequence[Operator],
    distance,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    simulator: BackendV2,
    shots_per_input: int,
):
    score_accumulator = 0
    sample_count = 0
    for init_operator in input_sampler:
        rho_a = to_bloch(zero.evolve(init_operator))
        measurements_xyz = run_simulation(
            init_operator, alice_noise, bob_noise, shots_per_input, simulator
        )
        p_is, rho_B_is = get_probabilities_and_states(measurements_xyz, shots_per_input)
        score_accumulator += get_score_protocol_distance(
            p_is, rho_a, rho_B_is, distance
        )
        sample_count += 1

    return score_accumulator / sample_count
