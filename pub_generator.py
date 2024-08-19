from typing import Iterable
from qiskit.primitives.containers.sampler_pub import SamplerPubLike
from qiskit.primitives.containers.bindings_array import BindingsArrayLike
from channels import BaseChannel
from teleportation_circuit import get_circuit
from math import pi
import numpy as np


to_z_basis_rotation = [0, 0, 0]
to_x_basis_rotation = [pi / 2, 0, pi]
to_y_basis_rotation = [pi / 2, 0, pi / 2]
basis_change = np.asarray(
    [to_x_basis_rotation, to_y_basis_rotation, to_z_basis_rotation]
)


def generate_pub(
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    exploration_space: Iterable[tuple[float, float]],
    input_sampler: np.ndarray,
    transpiler,
) -> SamplerPubLike:
    qc = transpiler(get_circuit(alice_noise, bob_noise))
    noise_params = np.array(
        [
            [
                pA := alice_noise.get_theta(x),
                pB := bob_noise.get_theta(y),
                BaseChannel.bob_optimal_rotation_for_noises(
                    alice_noise, bob_noise, pA, pB
                ),
            ]
            for x, y in exploration_space
        ]
    )

    cart_prod = get_cartesian_product_nashe(noise_params, input_sampler, basis_change)

    """ a = np.repeat(basis_change, len(input_sampler) * len(noise_params), axis=0)
    b = np.repeat(
        np.tile(input_sampler, (len(basis_change), 1)), len(noise_params), axis=0
    )
    c = np.tile(noise_params, (len(basis_change) * len(input_sampler), 1)) equiv, testear tiempos"""

    parameters: BindingsArrayLike = {
        ("theta_ADC_alice", "theta_ADC_bob", "bob_optimal_rotation_theta"): cart_prod[:, :3],
        ("init_theta", "init_phi", "init_lambda"): cart_prod[:, 3:6],
        ("basis_change_theta","basis_change_phi", "basis_change_lambda"): cart_prod[:, 6:],
    }
    return (transpiler(qc), parameters)


def get_cartesian_product_nashe(A, B, C):
    n = len(A)
    m = len(B)
    o = len(C)
    cart_prod = np.empty((n * m * o, 9), dtype='float64')
    for i in range(n):
        for j in range(m):
            for k in range(o):
                idx = i * m * o + j * o + k
                cart_prod[idx, :3] = A[i]
                cart_prod[idx, 3:6] = B[j]
                cart_prod[idx, 6:] = C[k]
    return cart_prod