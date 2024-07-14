from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit.library import UnitaryGate
def to_z_basis(qc, qubit):
    pass

def to_x_basis(qc, qubit):
    qc.h(qubit)

def to_y_basis(qc, qubit):
    qc.sdg(qubit)
    qc.h(qubit)
    
    
def run_simulation(initialization, alice_noise, bob_noise, shots, simulator):
    qc_x = get_circuit(initialization, alice_noise, bob_noise, to_x_basis)
    qc_y = get_circuit(initialization, alice_noise, bob_noise, to_y_basis)
    qc_z = get_circuit(initialization, alice_noise, bob_noise, to_z_basis)
    return run_simulation_in_basis(qc_x, shots, simulator), run_simulation_in_basis(qc_y, shots, simulator), run_simulation_in_basis(qc_z, shots, simulator)
    
def run_simulation_in_basis(qc, shots, simulator):
    circ = transpile(qc, simulator)
    # Run and get counts
    result = simulator.run(circ, shots=shots).result()
    return result.get_counts(circ)

def entangle(qc, q1, q2):
    qc.h(q1)
    qc.cx(q1, q2)
    qc.barrier()
    
def bell_basis(qc, q1, q2):
    qc.cx(q1, q2)
    qc.h(q1)
    qc.barrier()
    
def bob_correction(qc, q1, q2, qbob):
    qc.cx(q1, qbob)
    qc.cz(q2, qbob)
    
def get_circuit(init_operator, alice_noise, bob_noise, basis_change):
    alice_reg = QuantumRegister(2, name="alice")
    bob_reg = QuantumRegister(1, name="bob")
    qc = QuantumCircuit(alice_reg, bob_reg)
    bob_class_reg = ClassicalRegister(1)
    qc.add_register(bob_class_reg)
    alice_bell_class_reg = ClassicalRegister(1)
    qc.add_register(alice_bell_class_reg)
    alice_psi_class_reg = ClassicalRegister(1)
    qc.add_register(alice_psi_class_reg)
    
    qc.append(UnitaryGate(init_operator.data, label='State Initializator'), [alice_reg[0]])

    entangle(qc, alice_reg[1], bob_reg[0])
    alice_noise(qc, alice_reg[1])
    bob_noise(qc, bob_reg[0])
    bell_basis(qc, alice_reg[0], alice_reg[1])
    bob_correction(qc, alice_reg[1], alice_reg[0], bob_reg[0])

    basis_change(qc, bob_reg[0])
    
    qc.measure(bob_reg[0], bob_class_reg)
    qc.measure(alice_reg[0], alice_psi_class_reg)
    qc.measure(alice_reg[1], alice_bell_class_reg)
    return qc