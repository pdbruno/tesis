from typing import Sequence
from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit
from qiskit.circuit import Parameter
import numpy.polynomial.polynomial as poly
import numpy as np
import math
from abc import ABC, abstractmethod
from qiskit.circuit.operation import Operation
from qiskit.circuit.library import IGate, ZGate

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
        return self.get_circuit(suffix).to_instruction(label=self.get_instruction_label(suffix))

    def __init__(self) -> None:
        super().__init__()
        self.label = None

    @staticmethod
    def bob_optimal_rotation_for_noises(noise_a, noise_b, pA, pB) -> float:
        if (
            type(noise_a) == AmplitudeDampingChannel
            and type(noise_b) == MirroredAmplitudeDampingChannel
        ) or (
            type(noise_a) == MirroredAmplitudeDampingChannel
            and type(noise_b) == AmplitudeDampingChannel
        ):
            if (
                0 < pB
                and pB < 1
                and 1 + math.sqrt(1 + 2 * pB - 3 * pB**2) < 2 * pA + pB
            ) or (pB == 1 and pA > 0):
                """ 00 -> Z
                    01 -> ZX
                    10 -> I
                    11 -> X """
                return math.pi #ZGate

                """ 00 -> I
                    01 -> X
                    10 -> Z
                    11 -> ZX """
        return 0


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
        self.label = "DepC"

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
        qc = QuantumCircuit(2, 1)
        theta = Parameter(self.get_parameter_label(suffix))
        qc.cry(theta, 0, 1)
        qc.cx(1, 0)
        qc.measure(1, 0)
        return qc

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_zero = QuantumRegister(
            1, f"zero_{self.get_instruction_label(qubit.name)}"
        )
        parent_meas = ClassicalRegister(
            1, f"meas_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_zero)
        qc.add_register(parent_meas)
        qc.append(self.get_instruction(qubit.name), [qubit, parent_zero], [parent_meas])

    


class MirroredAmplitudeDampingChannel(AmplitudeDampingChannel):
    def get_circuit(self, suffix):
        qc = super().get_circuit(suffix)
        qc.x(0)
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
