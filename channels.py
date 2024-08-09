from typing import Sequence
from qiskit import QuantumRegister, QuantumCircuit
from qiskit.circuit import Parameter
import numpy.polynomial.polynomial as poly
import numpy as np
import math
from abc import ABC, abstractmethod
from qiskit.circuit.operation import Operation
from qiskit.circuit.library import IGate, ZGate

class BaseChannel(ABC):
    @staticmethod
    @abstractmethod
    def get_circuit() -> QuantumCircuit:
        pass

    @classmethod
    def get_instruction_label(cls, suffix) -> str:
        return f"{cls.label}_" + suffix

    @abstractmethod
    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        pass

    @abstractmethod
    def get_parameters(self) -> Sequence:
        pass

    def get_instruction(self, suffix):
        return self.circuit.assign_parameters(self.get_parameters()).to_instruction(
            label=self.get_instruction_label(suffix)
        )

    circuit = get_circuit()
    label = None

    @staticmethod
    def bob_optimal_rotation_for_noises(noise_a, noise_b) -> Operation:
        if (
            type(noise_a) == AmplitudeDampingChannel
            and type(noise_b) == MirroredAmplitudeDampingChannel
        ) or (
            type(noise_a) == MirroredAmplitudeDampingChannel
            and type(noise_b) == AmplitudeDampingChannel
        ):
            pA = noise_a.p
            pB = noise_b.p

            if (
                0 < pB
                and pB < 1
                and 1 + math.sqrt(1 + 2 * pB - 3 * pB**2) < 2 * pA + pB
            ) or (pB == 1 and pA > 0):
                """ 00 -> Z
                    01 -> ZX
                    10 -> I
                    11 -> X """
                return ZGate()

                """ 00 -> I
                    01 -> X
                    10 -> Z
                    11 -> ZX """
        return IGate()


class DepolarizingChannel(BaseChannel):

    def get_parameters(self) -> Sequence:
        return [self.get_theta()]

    @staticmethod
    def get_circuit():
        qc = QuantumCircuit(4)
        theta = Parameter("theta")
        qc.ry(theta, 1)
        qc.ry(theta, 2)
        qc.ry(theta, 3)
        qc.cx(1, 0)
        qc.cy(2, 0)
        qc.cz(3, 0)
        return qc

    def get_theta(self):
        z = poly.polyroots([-self.p / 4, 1, -1])[0]
        # p/4=z−z^2 --------- For any 0⩽p⩽1, there should be a unique 0⩽z⩽1/2 which solves the above equation.
        return 2 * math.acos(math.sqrt(1 - z))  # cos(theta/2) = sqrt(1-z)

    def __init__(self, p: float) -> None:
        super().__init__()
        self.p = p

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_delta_z = QuantumRegister(
            3, f"delta_z_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_delta_z)
        qc.append(
            self.get_instruction(qubit.name),
            [qubit, parent_delta_z[0], parent_delta_z[1], parent_delta_z[2]],
        )

    circuit = get_circuit()
    label = "DepC"


class AmplitudeDampingChannel(BaseChannel):

    def get_parameters(self) -> Sequence:
        return [self.get_theta()]

    def get_theta(self):
        return 2 * math.asin(math.sqrt(self.p))  # sin^2 (θ/2) = γ

    def __init__(self, p: float) -> None:
        super().__init__()
        self.p = p

    @staticmethod
    def get_circuit():
        qc = QuantumCircuit(2)
        theta = Parameter("theta")
        qc.cry(theta, 0, 1)
        qc.cx(1, 0)
        #qc.measure(1, 0)
        return qc

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        parent_zero = QuantumRegister(
            1, f"zero_{self.get_instruction_label(qubit.name)}"
        )
        qc.add_register(parent_zero)
        qc.append(self.get_instruction(qubit.name), [qubit, parent_zero])

    circuit = get_circuit()
    label = "ADC"


class MirroredAmplitudeDampingChannel(AmplitudeDampingChannel):

    @staticmethod
    def get_circuit():
        qc = AmplitudeDampingChannel.get_circuit()
        qc.x(0)
        return qc

    circuit = get_circuit()
    label = "MADC"


class NoiselessChannel(BaseChannel):
    @staticmethod
    def get_circuit() -> QuantumCircuit:  # type: ignore
        pass

    def append_channel_instruction(self, qc: QuantumCircuit, qubit: QuantumRegister):
        pass

    def get_parameters(self) -> Sequence:  # type: ignore
        pass
