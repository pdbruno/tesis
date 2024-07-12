from qiskit import QuantumRegister
import numpy.polynomial.polynomial as poly
import math

def depolarizing_channel(p):
    def get_angle_for_probability(p):
        z = poly.polyroots([-p/4, 1, -1])[0] #p/4=z−z^2 --------- For any 0⩽p⩽1, there should be a unique 0⩽z⩽1/2 which solves the above equation.
        return 2*math.acos(math.sqrt(1-z)) #cos(theta/2) = sqrt(1-z)
    
    def channel(qc, qubit):
        delta_z = QuantumRegister(3)
        qc.add_register(delta_z)
        angle = get_angle_for_probability(p)
        qc.ry(angle, delta_z[0])
        qc.ry(angle, delta_z[1])
        qc.ry(angle, delta_z[2])
        qc.cx(delta_z[0], qubit)
        qc.cy(delta_z[1], qubit)
        qc.cz(delta_z[2], qubit)
    return channel
    
def amplitude_damping_channel(p):
    def get_angle_for_probability(p):
        return math.asin(math.sqrt(p)) #γ = sin^2(θ)
    
    def channel(qc, qubit):
        zero = QuantumRegister(1)
        qc.add_register(zero)
        angle = get_angle_for_probability(p)
        qc.cry(angle, qubit, zero[0])
        qc.cx(zero[0], qubit)
    return channel

def mirrored_amplitude_damping_channel(p):
    def channel(qc, qubit):
        amplitude_damping_channel(p)(qc, qubit)
        qc.x(qubit)
    return channel