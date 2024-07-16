from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator
from channels import BaseChannel
from qiskit.providers import BackendV2

x_basis_change = Operator.from_label("H")
y_basis_change = Operator.from_label("S").adjoint() & Operator.from_label("H")
z_basis_change = Operator.from_label("I")


def to_z_basis(qc, qubit):
    pass


def to_x_basis(qc, qubit):
    qc.h(qubit)


def to_y_basis(qc, qubit):
    qc.sdg(qubit)
    qc.h(qubit)


def run_simulation(
    init_operator: Operator,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    shots: int,
    simulator: BackendV2,
):
    qc_x = get_circuit(init_operator, alice_noise, bob_noise, x_basis_change)
    qc_y = get_circuit(init_operator, alice_noise, bob_noise, y_basis_change)
    qc_z = get_circuit(init_operator, alice_noise, bob_noise, z_basis_change)
    return (
        run_simulation_in_basis(qc_x, shots, simulator),
        run_simulation_in_basis(qc_y, shots, simulator),
        run_simulation_in_basis(qc_z, shots, simulator),
    )


def run_simulation_in_basis(qc: QuantumCircuit, shots: int, simulator: BackendV2):
    circ = transpile(qc, simulator)
    # Run and get counts
    result = simulator.run(circ, shots=shots).result()  # type: ignore
    return result.get_counts(circ)


def entangle(qc, q1, q2):
    qc.h(q1)
    qc.cx(q1, q2)


def bell_basis(qc, q1, q2):
    qc.cx(q1, q2)
    qc.h(q1)


def bob_correction(
    qc: QuantumCircuit, q1, q2, qbob, alice_noise: BaseChannel, bob_noise: BaseChannel
):
    # aca va otra correccion que depende de los dos ruidos
    qc.append(
        UnitaryGate(BaseChannel.bob_optimal_rotation_for_noises(alice_noise, bob_noise), label="Optimal Rotation"), [qbob]  # type: ignore
    )
    qc.cx(q1, qbob)
    qc.cz(q2, qbob)


def get_circuit(
    init_operator: Operator,
    alice_noise: BaseChannel,
    bob_noise: BaseChannel,
    basis_change: Operator,
):
    input_reg = QuantumRegister(1, name="input")
    alice_reg = QuantumRegister(1, name="alice")
    bob_reg = QuantumRegister(1, name="bob")
    qc = QuantumCircuit(input_reg, alice_reg, bob_reg)
    bob_class_reg = ClassicalRegister(1)
    qc.add_register(bob_class_reg)
    alice_bell_class_reg = ClassicalRegister(1)
    qc.add_register(alice_bell_class_reg)
    alice_psi_class_reg = ClassicalRegister(1)
    qc.add_register(alice_psi_class_reg)

    qc.append(
        UnitaryGate(init_operator.data, label="State Initializator"), [input_reg[0]]  # type: ignore
    )

    entangle(qc, alice_reg[0], bob_reg[0])
    qc.barrier()
    alice_noise.append_channel_instruction(qc, alice_reg)
    bob_noise.append_channel_instruction(qc, bob_reg)
    qc.barrier()
    bell_basis(qc, input_reg[0], alice_reg[0])
    qc.barrier()
    bob_correction(qc, alice_reg[0], input_reg[0], bob_reg[0], alice_noise, bob_noise)

    qc.barrier()
    qc.append(
        UnitaryGate(basis_change.data, label="Basis Change"), [bob_reg[0]]  # type: ignore
    )
    qc.barrier()

    qc.measure(bob_reg[0], bob_class_reg)
    qc.measure(input_reg[0], alice_psi_class_reg)
    qc.measure(alice_reg[0], alice_bell_class_reg)
    return qc
