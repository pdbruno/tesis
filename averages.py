import numpy as np
from quantum_state_tomography import get_probabilities_and_states, quantum_state_tomography
from math import cos, sin
from cmath import exp

""" 
           ( cos(θ/2)          -e^(iλ)*sin(θ/2) )
U(θ,ϕ,λ) = (                                    )
           ( e^(iϕ)*sin(θ/2)   e^(ϕ+λ)*cos(θ/2) )

 """
def to_bloch(eulerian_angles):
    [theta, phi, lam] = eulerian_angles
    return np.array(
        [
            0.5 * (exp(-phi * 1j) * (1 + exp(phi * 2j)) * sin(theta)).real,
            0.5 * (exp(-phi * 1j) * (-1 + exp(phi * 2j)) * sin(theta)).imag,
            cos(theta),
        ]
    )



def get_score_protocol_distance(p_is, rho_a, rho_B_is, distance):
    return p_is @ distance(rho_a, rho_B_is)



def average_distance_of_teleportation(result_for_input_state, rho_A, shots, distances):
    p_is, rho_B_is = quantum_state_tomography(result_for_input_state, shots)
    bloch_rho_A = to_bloch(rho_A)
    return [get_score_protocol_distance(p_is, bloch_rho_A, rho_B_is, d) for d in distances]

