from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit.library import UGate
from qiskit.circuit import Parameter, ParameterVector
from channels import BaseChannel


def entangle(qc, q1, q2):
    qc.h(q1)
    qc.cx(q1, q2)


def bell_basis(qc, q1, q2):
    qc.cx(q1, q2)
    qc.h(q1)


def bob_correction(qc: QuantumCircuit, q1, q2, qbob):
    qc.cx(q1, qbob)
    qc.cz(q2, qbob)
    # aca va otra correccion que depende de los dos ruidos
    qc.append(
        UGate(Parameter("bob_optimal_rotation_theta"), Parameter("bob_optimal_rotation_phi"), Parameter("bob_optimal_rotation_lambda"), label="Optimal Rotation"),
        [qbob],
    )


def get_circuit(alice_noise: BaseChannel, bob_noise: BaseChannel):
    input_reg = QuantumRegister(1, name="input")
    alice_reg = QuantumRegister(1, name="alice")
    bob_reg = QuantumRegister(1, name="bob")
    input_meas = ClassicalRegister(1, name="input_meas")
    alice_meas = ClassicalRegister(1, name="alice_meas")
    bob_meas = ClassicalRegister(1, name="bob_meas")

    qc = QuantumCircuit(input_reg, alice_reg, bob_reg, input_meas, alice_meas, bob_meas)

    initialize_input_state(input_reg, qc)

    entangle(qc, alice_reg[0], bob_reg[0])

    alice_noise.append_channel_instruction(qc, alice_reg)
    bob_noise.append_channel_instruction(qc, bob_reg)

    bell_basis(qc, input_reg[0], alice_reg[0])

    bob_correction(qc, alice_reg[0], input_reg[0], bob_reg[0])

    bob_basis_change(bob_reg, qc)

    qc.measure(bob_reg[0], bob_meas)
    qc.measure(input_reg[0], input_meas)
    qc.measure(alice_reg[0], alice_meas)
    return qc


def bob_basis_change(bob_reg, qc):
    qc.append(
        UGate(
            Parameter("basis_change_theta"),
            Parameter("basis_change_phi"),
            Parameter("basis_change_lambda"),
            label="Basis Change",
        ),
        [bob_reg[0]],
    )


def initialize_input_state(input_reg, qc):
    qc.append(
        UGate(
            Parameter("init_theta"),
            Parameter("init_phi"),
            Parameter("init_lambda"),
            label="State Initializator",
        ),
        [input_reg[0]],
    )