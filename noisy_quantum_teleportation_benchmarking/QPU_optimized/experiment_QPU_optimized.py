from noisy_quantum_teleportation_benchmarking.averages import get_teleportation_score_for_state
from noisy_quantum_teleportation_benchmarking.channels import BaseChannel
import numpy as np
from qiskit import QuantumCircuit
from typing import Callable, Sequence
from noisy_quantum_teleportation_benchmarking.distances import DistanceFunction
from noisy_quantum_teleportation_benchmarking.sampler import PauliSampler
from noisy_quantum_teleportation_benchmarking.QPU_optimized.teleportation_circuit_QPU_optimized import get_circuits
from qiskit.primitives.base import BaseSamplerV2

class ExperimentQPUOptimized:
    def __init__(
        self,
        channel_combinations: Sequence[tuple[BaseChannel, BaseChannel]],
        exploration_space: list[tuple[float, float]],
        transpiler: Callable[[QuantumCircuit], QuantumCircuit],
    ) -> None:
        self.channel_combinations = channel_combinations
        self.exploration_space = exploration_space
        self.transpiler = transpiler
        self.pubs = [
            pub
            for chA, chB in channel_combinations
            for pub in self.generate_pubs(chA, chB)
        ]

    def generate_pubs(
        self, alice_noise: BaseChannel, bob_noise: BaseChannel
    ) -> list[QuantumCircuit]:
        return [
            self.transpiler(qc)
            for pA, pB in self.exploration_space
            for qc in get_circuits(alice_noise, bob_noise, pA, pB)
        ]

    def run_with_sampler(
        self,
        sampler: BaseSamplerV2,
        distances: list[DistanceFunction],
        shots: int,
        file_path: str,
    ):
        results = sampler.run(self.pubs, shots=shots).result()
        with open(file_path, "w+") as f:
            f.write("chA,chB,pA,pB,d,score\n")
            i = 0
            for chA, chB in self.channel_combinations:
                for pA, pB in self.exploration_space:
                    scores = []
                    for phi_eulerian_angles in PauliSampler.all:
                        result_for_input_state = []
                        for meas_basis in range(3):
                            pub_result = results[i]
                            i += 1
                            measurements = np.array(
                                [
                                    np.squeeze(pub_result.data["input_meas"].array),
                                    np.squeeze(pub_result.data["alice_meas"].array),
                                    np.squeeze(pub_result.data["bob_meas"].array),
                                ]
                            )
                            result_for_input_state.append(measurements)

                        scores.append(
                            get_teleportation_score_for_state(
                                result_for_input_state,
                                phi_eulerian_angles,
                                shots,
                                distances,
                            )
                        )

                    avg_distances = np.array(scores).mean(axis=0)
                    f.writelines(
                        f"{chA.label},{chB.label},{pA},{pB},{d.__name__},{avg_distance_of_teleportation}\n"
                        for d, avg_distance_of_teleportation in zip(
                            distances, avg_distances
                        )
                    )
