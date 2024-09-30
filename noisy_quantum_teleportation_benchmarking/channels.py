from __future__ import annotations
from typing import Any
from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit
from qiskit.circuit import Parameter
import numpy.polynomial.polynomial as poly
import numpy as np
import math
from abc import ABC, abstractmethod

bell_states = [
    [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
    [0, 1 / math.sqrt(2), 1 / math.sqrt(2), 0],
    [1 / math.sqrt(2), 0, 0, -(1 / math.sqrt(2))],
    [0, 1 / math.sqrt(2), -(1 / math.sqrt(2)), 0],
]

# 0: I, 1: X, 2: Z, 3: ZX

optimal_rotation_angles = [
    [0, 0, 0],
    [np.pi, 0, np.pi],
    [0, 0, np.pi],
    [np.pi, np.pi, np.pi],
]


class BaseChannel(ABC):
    @abstractmethod
    def get_circuit(self, suffix) -> QuantumCircuit:
        pass

    @abstractmethod
    def get_theta(self, p) -> float:
        pass

    def get_instruction_label(self, suffix) -> str:
        return f"{self.label}_{suffix}"

    def get_parameter_label(self, suffix) -> str:
        return f"theta_{self.label}_{suffix}"

    @abstractmethod
    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        pass

    def get_instruction(self, suffix):
        return self.get_circuit(suffix).to_instruction(
            label=self.get_instruction_label(suffix)
        )

    def __init__(self) -> None:
        super().__init__()
        self.label = ""

    @staticmethod
    def l_max_for_noises(
        noise_a: BaseChannel, noise_b: BaseChannel, pA: float, pB: float
    ) -> int:

        eigenvalues, eigenvectors = eigencosas[(noise_a.label, noise_b.label)](pA, pB)
        max_eigenvalue = max(eigenvalues)
        max_eigenvectors = [
            eigenvectors[i]
            for i, eigenvalue in enumerate(eigenvalues)
            if eigenvalue == max_eigenvalue
        ]
        matching_bell_state = next(
            (x for x in max_eigenvectors if x in bell_states), None
        )  # first that matches
        if matching_bell_state == None:
            matching_bell_state = np.sum(
                max_eigenvectors, axis=0
            ).tolist()  # linear combination
            if matching_bell_state not in bell_states:
                matching_bell_state = bell_states[0]  # anyone

        return bell_states.index(matching_bell_state) + 1

    @staticmethod
    def bob_optimal_rotation_for_noises(
        noise_a: BaseChannel, noise_b: BaseChannel, pA: float, pB: float
    ) -> list[float]:
        return optimal_rotation_angles[
            BaseChannel.l_max_for_noises(noise_a, noise_b, pA, pB) - 1
        ]


class DepolarizingChannel(BaseChannel):

    def get_circuit(self, suffix):
        qc = QuantumCircuit(4)
        theta = Parameter(self.get_parameter_label(suffix))
        qc.ry(theta, 1)
        qc.ry(theta, 2)
        qc.ry(theta, 3)
        qc.cx(1, 0)
        qc.cy(2, 0)
        qc.cz(3, 0)
        return qc

    def get_theta(self, p):
        z = poly.polyroots([-p / 4, 1, -1])[0]
        # p/4=z−z^2 --------- For any 0⩽p⩽1, there should be a unique 0⩽z⩽1/2 which solves the above equation.
        return 2 * math.acos(math.sqrt(1 - z))  # cos(theta/2) = sqrt(1-z)

    def __init__(self) -> None:
        super().__init__()
        self.label = "DC"

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_delta_z = QuantumRegister(
            3, f"delta_z_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_delta_z)
        qc.append(
            self.get_instruction(qubit.name),
            [qubit, parent_delta_z[0], parent_delta_z[1], parent_delta_z[2]],
        )


class AmplitudeDampingChannel(BaseChannel):

    def get_theta(self, p):
        return 2 * math.asin(math.sqrt(p))  # sin^2 (θ/2) = γ

    def __init__(self) -> None:
        super().__init__()
        self.label = "ADC"

    def get_circuit(self, suffix):
        qc = QuantumCircuit(2)
        theta = Parameter(self.get_parameter_label(suffix))
        qc.cry(theta, 0, 1)
        qc.cx(1, 0)
        return qc

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_zero = QuantumRegister(
            1, f"zero_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_zero)
        qc.append(self.get_instruction(qubit.name), [qubit, parent_zero])


class PhaseDampingChannel(BaseChannel):

    def get_theta(self, p):
        return 2 * math.asin(math.sqrt(p))  # sin^2 (θ/2) = γ
        return 2 * math.asin(p)  # sin (θ/2) = γ

    def __init__(self) -> None:
        super().__init__()
        self.label = "PDC"

    def get_circuit(self, suffix):
        qc = QuantumCircuit(2)
        theta = Parameter(self.get_parameter_label(suffix))
        qc.cry(theta, 0, 1)
        return qc

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_zero = QuantumRegister(
            1, f"zero_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_zero)
        qc.append(self.get_instruction(qubit.name), [qubit, parent_zero])


class MirroredAmplitudeDampingChannel(AmplitudeDampingChannel):
    def get_circuit(self, suffix):
        qc = QuantumCircuit(2)
        theta = Parameter(self.get_parameter_label(suffix))
        qc.x(0)
        qc.x(1)
        qc.cry(theta, 0, 1)
        qc.cx(1, 0)
        return qc

    def __init__(self) -> None:
        super().__init__()
        self.label = "MADC"


class NoiselessChannel(BaseChannel):
    def get_circuit(self, suffix):
        pass

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        pass

    def get_theta(self, p) -> float:
        return 0

    def __init__(self) -> None:
        super().__init__()
        self.label = None


eigencosas = {
    ("ADC", "ADC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -2 * pA * pB),
            0.25 * (pA + pB + -2 * pA * pB),
            0.25
            * (2 + -1 * pB + pA * (-1 + 2 * pB) + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25
            * (2 + -1 * pB + pA * (-1 + 2 * pB) + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("ADC", "MADC"): lambda pA, pB: [
        [
            0.25 * (2 + -1 * pA + -1 * pB + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pA + -1 * pB + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (pA + pB),
            0.25 * (pA + pB),
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
        ],
    ],
    ("ADC", "DC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + 2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)
            ),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + -2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)
            ),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("ADC", "PDC"): lambda pA, pB: [
        [
            0.25 * pA,
            0.25 * pA,
            0.25 * (2 + -1 * pA + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pA + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("MADC", "ADC"): lambda pA, pB: [
        [
            0.25 * (2 + -1 * pA + -1 * pB + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pA + -1 * pB + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (pA + pB),
            0.25 * (pA + pB),
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
        ],
    ],
    ("MADC", "MADC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -2 * pA * pB),
            0.25 * (pA + pB + -2 * pA * pB),
            0.25
            * (2 + -1 * pB + pA * (-1 + 2 * pB) + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25
            * (2 + -1 * pB + pA * (-1 + 2 * pB) + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("MADC", "DC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + 2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)
            ),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + -2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)
            ),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("MADC", "PDC"): lambda pA, pB: [
        [
            0.25 * pA,
            0.25 * pA,
            0.25 * (2 + -1 * pA + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pA + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("DC", "ADC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + 2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)
            ),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + -2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)
            ),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("DC", "MADC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + 2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)
            ),
            0.25
            * (
                2
                + -1 * pA
                + -1 * pB
                + pA * pB
                + -2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)
            ),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("DC", "DC"): lambda pA, pB: [
        [
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (pA + pB + -1 * pA * pB),
            0.25 * (4 + -3 * pB + 3 * pA * (-1 + pB)),
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("DC", "PDC"): lambda pA, pB: [
        [
            0.25 * pA,
            0.25 * pA,
            0.25 * (2 + -1 * pA + 2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)),
            0.25 * (2 + -1 * pA + -2 * (-1 + pA) * pow((1 + -1 * pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
    ("PDC", "ADC"): lambda pA, pB: [
        [
            0.25 * (2 + -1 * pB + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pB + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * pB,
            0.25 * pB,
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
        ],
    ],
    ("PDC", "MADC"): lambda pA, pB: [
        [
            0.25 * (2 + -1 * pB + -2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * (2 + -1 * pB + 2 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.25 * pB,
            0.25 * pB,
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
        ],
    ],
    ("PDC", "DC"): lambda pA, pB: [
        [
            0.25 * (2 + -1 * pB + 2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)),
            0.25 * (2 + -1 * pB + -2 * (-1 + pB) * pow((1 + -1 * pA), 0.5)),
            0.25 * pB,
            0.25 * pB,
        ],
        [
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
        ],
    ],
    ("PDC", "PDC"): lambda pA, pB: [
        [
            0,
            0,
            0.5 * (1 + -1 * pow((-1 + pA) * (-1 + pB), 0.5)),
            0.5 * (1 + pow((-1 + pA) * (-1 + pB), 0.5)),
        ],
        [
            [0, 0, 1 / math.sqrt(2), 0],
            [0, 1 / math.sqrt(2), 0, 0],
            [-1 * 1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
            [1 / math.sqrt(2), 0, 0, 1 / math.sqrt(2)],
        ],
    ],
}
