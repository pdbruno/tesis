from averages import average_distance_of_teleportation
from channels import BaseChannel
import numpy as np
from numpy.typing import NDArray
from typing import Callable, Tuple
from qiskit import QuantumCircuit
from typing import cast, Callable, Sequence
from qiskit.primitives.containers.sampler_pub import SamplerPubLike
from qiskit.primitives.containers.bindings_array import BindingsArrayLike
from channels import BaseChannel
from distances import DistanceFunction
from sampler import BaseInputSampler
from teleportation_circuit import (
    BASIS_CHANGE,
    BOB_OPTIMAL_ROTATION,
    INIT_STATE,
    get_circuit,
)
from math import pi
from qiskit.primitives.base import BaseSamplerV2

to_z_basis_rotation = [0, 0, 0]
to_x_basis_rotation = [pi / 2, 0, pi]
to_y_basis_rotation = [pi / 2, 0, pi / 2]
basis_change = np.asarray(
    [to_x_basis_rotation, to_y_basis_rotation, to_z_basis_rotation]
)


def get_cartesian_product_nashe(
    A: NDArray[np.floating], sampler: BaseInputSampler, C: NDArray[np.floating]
):
    n = len(A)
    m = sampler.length
    o = len(C)
    cart_prod = np.empty((n * m * o, 11), dtype="float64")
    for i in range(n):
        B = sampler.get_samples()
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
        channel_combinations: Sequence[tuple[BaseChannel, BaseChannel]],
        exploration_space: list[tuple[float, float]],
        input_sampler: BaseInputSampler,
        transpiler: Callable[[QuantumCircuit], QuantumCircuit],
    ) -> None:
        self.channel_combinations = channel_combinations
        self.exploration_space = exploration_space
        self.input_sampler = input_sampler
        self.transpiler = transpiler
        self.pubs = (self.generate_pub(chA, chB) for chA, chB in channel_combinations)

    def generate_pub(
        self, alice_noise: BaseChannel, bob_noise: BaseChannel
    ) -> Tuple[QuantumCircuit, BindingsArrayLike]:
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
                alice_noise.get_parameter_label("alice"),
                bob_noise.get_parameter_label("bob"),
                *BOB_OPTIMAL_ROTATION,
            ): cart_prod[:, :5],
            tuple(INIT_STATE): cart_prod[:, 5:8],
            tuple(BASIS_CHANGE): cart_prod[:, 8:],
        }
        return (self.transpiler(get_circuit(alice_noise, bob_noise)), parameters)

    def run_with_sampler(
        self,
        sampler: BaseSamplerV2,
        distances: list[DistanceFunction],
        shots: int,
        file_path: str,
    ):

        with open(file_path, "w+") as f:
            f.write("chA,chB,pA,pB,d,score\n")
            for (chA, chB), pub in zip(self.channel_combinations, self.pubs):
                pub_result = sampler.run([pub]).result()[0]
                for result_for_noise_configuration, (pA, pB), input_samples in zip(
                    self.reshape_result_data(pub_result),
                    self.exploration_space,
                    self.get_samples_for_noise_conf(pub[1]),
                ):
                    scores = self.get_scores_for_distances(
                        distances, shots, result_for_noise_configuration, input_samples
                    )

                    f.writelines(
                        f"{chA.label},{chB.label},{pA},{pB},{d.__name__},{score}\n"
                        for d, score in zip(distances, scores)
                    )

    def get_samples_for_noise_conf(self, bindings_array):
        params = bindings_array[tuple(INIT_STATE)]
        return np.split(
            params[::3, :], cast(int, len(params) / 3 / self.input_sampler.length)
        )

    def get_scores_for_distances(
        self, distances, shots, result_for_noise_configuration, input_samples
    ) -> NDArray[np.floating]:
        return np.array(
            [
                average_distance_of_teleportation(
                    result_for_input_state, phi_eulerian_angles, shots, distances
                )
                for result_for_input_state, phi_eulerian_angles in zip(
                    result_for_noise_configuration, input_samples
                )
            ]
        ).mean(axis=0)

    def reshape_result_data(self, pub_result):
        length = pub_result.data["input_meas"].shape[0]
        measurements = np.array(
            [
                np.squeeze(pub_result.data["input_meas"].array),
                np.squeeze(pub_result.data["alice_meas"].array),
                np.squeeze(pub_result.data["bob_meas"].array),
            ]
        )  # (3=measured_qubits, length, shots)

        measurements = np.swapaxes(
            measurements, 0, 1
        )  # (length, 3=measured_qubits, shots)

        by_samples = np.array(
            np.split(measurements, cast(int, length / 3))
        )  # (length / 3=measurement_basis, 3=measurement_basis, 3=measured_qubits, shots)

        by_noise_configuration = np.split(
            by_samples, cast(int, len(by_samples) / self.input_sampler.length)
        )
        # (121, input_samples, 3, 3, shots) = (noise_grid, input_samples, measurement_basis, measured_qubits, shots)
        return by_noise_configuration
