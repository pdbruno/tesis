from typing import Callable
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from noisy_quantum_teleportation_benchmarking.channels import BaseChannel


def to_z_basis(qc, qubit):
    pass


def to_x_basis(qc, qubit):
    qc.h(qubit)


def to_y_basis(qc, qubit):
    qc.sdg(qubit)
    qc.h(qubit)


def initialize_eigenvector_x_mas(qc, qubit):
    qc.h(qubit)


def initialize_eigenvector_x_menos(qc, qubit):
    qc.x(qubit)
    qc.h(qubit)


def initialize_eigenvector_y_mas(qc, qubit):
    qc.h(qubit)
    qc.s(qubit)


def initialize_eigenvector_y_menos(qc, qubit):
    qc.h(qubit)
    qc.sdg(qubit)


def initialize_eigenvector_z_mas(qc, qubit):
    pass


def initialize_eigenvector_z_menos(qc, qubit):
    qc.x(qubit)


type GateApplication = Callable[[QuantumCircuit, QuantumRegister], None]

basis: list[GateApplication] = [to_x_basis, to_y_basis, to_z_basis]

initializers: list[GateApplication] = [
    initialize_eigenvector_x_mas,
    initialize_eigenvector_x_menos,
    initialize_eigenvector_y_mas,
    initialize_eigenvector_y_menos,
    initialize_eigenvector_z_mas,
    initialize_eigenvector_z_menos,
]


def get_circuits(
    alice_noise: BaseChannel, bob_noise: BaseChannel, pA: float, pB: float
):
    return [
        get_circuit(init, alice_noise, bob_noise, basis_change, pA, pB)
        for init in initializers
        for basis_change in basis
    ]


def entangle(qc, q1, q2):
    qc.h(q1)
    qc.cx(q1, q2)


def bell_basis(qc, q1, q2):
    qc.cx(q1, q2)
    qc.h(q1)


def bob_correction(
    qc: QuantumCircuit,
    q1,
    q2,
    qbob,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    pA: float,
    pB: float,
):
    qc.cx(q1, qbob)
    qc.cz(q2, qbob)
    # aca va otra correccion que depende de los dos ruidos
    l_max = BaseChannel.l_max_for_noises(alice_noise, bob_noise, pA, pB)

    match l_max:
        case 1:
            pass
        case 2:
            qc.x(qbob)
        case 3:
            qc.z(qbob)
        case 4:
            qc.x(qbob)
            qc.z(qbob)


def get_circuit(
    initializer: GateApplication,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    basis_change: GateApplication,
    pA: float,
    pB: float,
) -> QuantumCircuit:
    input_reg = QuantumRegister(1, name="input")
    alice_reg = QuantumRegister(1, name="alice")
    bob_reg = QuantumRegister(1, name="bob")
    input_meas = ClassicalRegister(1, name="input_meas")
    alice_meas = ClassicalRegister(1, name="alice_meas")
    bob_meas = ClassicalRegister(1, name="bob_meas")

    qc = QuantumCircuit(input_reg, alice_reg, bob_reg, input_meas, alice_meas, bob_meas)

    initializer(qc, input_reg)

    entangle(qc, alice_reg[0], bob_reg[0])

    alice_noise.append_channel_instruction(qc, alice_reg)
    bob_noise.append_channel_instruction(qc, bob_reg)

    bell_basis(qc, input_reg[0], alice_reg[0])

    bob_correction(
        qc, alice_reg[0], input_reg[0], bob_reg[0], alice_noise, bob_noise, pA, pB
    )

    basis_change(qc, bob_reg)

    qc.measure(bob_reg[0], bob_meas)
    qc.measure(input_reg[0], input_meas)
    qc.measure(alice_reg[0], alice_meas)
    return qc.assign_parameters(
        {
            alice_noise.get_parameter_label("alice"): alice_noise.get_theta(pA),
            bob_noise.get_parameter_label("bob"): bob_noise.get_theta(pB),
        },  # type: ignore
        flat_input=False,
        inplace=False
    )
