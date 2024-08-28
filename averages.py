from typing import Iterable, List, Tuple
import random
import numpy as np
from channels import BaseChannel
from teleportation_circuit import get_circuit
from quantum_state_tomography import get_probabilities_and_states
from qiskit.providers import BackendV2
from math import cos, sin, pi
from cmath import exp
from numpy.typing import ArrayLike
from scipy.stats import rv_continuous

class sin_prob_dist(rv_continuous):
    def _pdf(self, theta):
        # The 0.5 is so that the distribution is normalized
        return 0.5 * np.sin(theta)

# Samples of theta should be drawn from between 0 and pi
sin_sampler = sin_prob_dist(a=0, b=np.pi)


def to_bloch(eulerian_angles):
    [theta, phi, lam] = eulerian_angles
    return np.array(
        [
            0.5 * (exp(-phi * 1j) * (1 + exp(phi * 2j)) * sin(theta)).real,
            0.5 * (exp(-phi * 1j) * (-1 + exp(phi * 2j)) * sin(theta)).imag,
            cos(theta),
        ]
    )


def haar_measure_average(shots) -> np.ndarray:
    thetas = sin_sampler.rvs(size=shots)
    phis_lambdas = np.random.uniform(0, 2 * pi, (shots, 2))
    raise ValueError()
    return np.concatenate([thetas, phis_lambdas], axis=1)


""" 
           ( cos(θ/2)          -e^(iλ)*sin(θ/2) )
U(θ,ϕ,λ) = (                                    )
           ( e^(iϕ)*sin(θ/2)   e^(ϕ+λ)*cos(θ/2) )

 """
eigenvector_x_mas = [pi / 2, 0, pi]  # [1., 0., 0.]
eigenvector_x_menos = [-pi / 2, 0, pi]  # [-1., 0., 0.]
eigenvector_y_mas = [-pi / 2, -pi / 2, pi / 2]  # [0., 1., 0.]
eigenvector_y_menos = [pi / 2, -pi / 2, pi / 2]  # [0., -1., 0.]
eigenvector_z_mas = [0.0, 0.0, 0.0]  # [0., 0., 1.]
eigenvector_z_menos = [pi, 0, 0]  # [0., 0., -1.]

all = np.array(
    [
        eigenvector_x_mas,
        eigenvector_x_menos,
        eigenvector_y_mas,
        eigenvector_y_menos,
        eigenvector_z_mas,
        eigenvector_z_menos,
    ]
)


def pauli_eigenvectors_average() -> np.ndarray:
    return all


def get_score_protocol_distance(p_is, rho_a, rho_B_is, distance):
    return p_is @ distance(rho_a, rho_B_is)


""" def average_distance_of_teleportation(
    input_sampler: Iterable[Tuple[float, float, float]],
    distance,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    simulator: BackendV2,
    shots_per_input: int,
):
    score_accumulator = 0
    sample_count = 0
    for theta, phi, lam in input_sampler:
        rho_a = to_bloch(theta, phi, lam)
        circuit = get_circuit(alice_noise, bob_noise)

        measurements_xyz = run_simulation(
            eigenvector_y_mas, alice_noise, bob_noise, shots_per_input, simulator
        )
        p_is, rho_B_is = get_probabilities_and_states(measurements_xyz, shots_per_input)
        score_accumulator += get_score_protocol_distance(
            p_is, rho_a, rho_B_is, distance
        )
        sample_count += 1

    return score_accumulator / sample_count """
