from averages import average_distance_of_teleportation
from channels import BaseChannel
import numpy as np
from numpy.typing import NDArray
from typing import Callable
from qiskit import QuantumCircuit
from typing import cast
from qiskit.primitives.containers.sampler_pub import SamplerPubLike
from qiskit.primitives.containers.bindings_array import BindingsArrayLike
from channels import BaseChannel
from teleportation_circuit import get_circuit
from math import pi
from qiskit.primitives.base import BaseSamplerV2
import pandas as pd

to_z_basis_rotation = [0, 0, 0]
to_x_basis_rotation = [pi / 2, 0, pi]
to_y_basis_rotation = [pi / 2, 0, pi / 2]
basis_change = np.asarray(
    [to_x_basis_rotation, to_y_basis_rotation, to_z_basis_rotation]
)


def get_cartesian_product_nashe(A, B, C):
    n = len(A)
    m = len(B)
    o = len(C)
    cart_prod = np.empty((n * m * o, 11), dtype="float64")
    for i in range(n):
        for j in range(m):
            for k in range(o):
                idx = i * m * o + j * o + k
                cart_prod[idx, :5] = A[i]
                cart_prod[idx, 5:8] = B[j]
                cart_prod[idx, 8:] = C[k]
    return cart_prod


class Experiment:
    def __init__(
        self,
        channel_combinations: list[tuple[BaseChannel, BaseChannel]],
        exploration_space: list[tuple[float, float]],
        input_sampler: NDArray[np.floating],
        transpiler: Callable[[QuantumCircuit], QuantumCircuit],
    ) -> None:
        self.channel_combinations = channel_combinations
        self.exploration_space = exploration_space
        self.input_sampler = input_sampler
        self.transpiler = transpiler
        self.pubs = [
            self.generate_pub(chA, chB)
            for chA, chB in channel_combinations
        ]

    def generate_pub(
        self, alice_noise: BaseChannel, bob_noise: BaseChannel
    ) -> SamplerPubLike:
        noise_params = np.array(
            [
                [
                    alice_noise.get_theta(pA),
                    bob_noise.get_theta(pB),
                    *BaseChannel.bob_optimal_rotation_for_noises(
                        alice_noise, bob_noise, pA, pB
                    ),
                ]
                for pA, pB in self.exploration_space
            ]
        )

        cart_prod = get_cartesian_product_nashe(
            noise_params, self.input_sampler, basis_change
        )

        """ a = np.repeat(basis_change, len(input_sampler) * len(noise_params), axis=0)
        b = np.repeat(
            np.tile(input_sampler, (len(basis_change), 1)), len(noise_params), axis=0
        )
        c = np.tile(noise_params, (len(basis_change) * len(input_sampler), 1)) equiv, testear tiempos"""

        parameters: BindingsArrayLike = {
            (
                "theta_ADC_alice",
                "theta_ADC_bob",
                "bob_optimal_rotation_theta",
                "bob_optimal_rotation_phi",
                "bob_optimal_rotation_lambda",
            ): cart_prod[:, :5],
            ("init_theta", "init_phi", "init_lambda"): cart_prod[:, 5:8],
            (
                "basis_change_theta",
                "basis_change_phi",
                "basis_change_lambda",
            ): cart_prod[:, 8:],
        }
        return (self.transpiler(get_circuit(alice_noise, bob_noise)), parameters)

    def run_with_sampler(self, sampler: BaseSamplerV2, distances, shots):
        results = sampler.run(self.pubs).result()

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

        for pub_result, (chA, chB) in zip(results, self.channel_combinations):
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
            # el 6 deberia reemplazarlo por la cantidad de input samples
            # (121, 6, 3, 3, 1000) = (noise_grid, input_sample, measurement_basis, measured_qubits, shots)


            for res_for_noise_configuration, (pA, pB) in zip(
                formatted_results, self.exploration_space
            ):
                scores = np.array(
                    [
                        average_distance_of_teleportation(
                            res_for_input_state, phi_eulerian_angles, shots, distances
                        )
                        for res_for_input_state, phi_eulerian_angles in zip(
                            res_for_noise_configuration, self.input_sampler
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

        return df