from numpy import linalg as LA
import math

def trace_distance(r, s):
    return LA.vector_norm(r-s);

def fidelity(r, s):
    return (1 + r @ s)/2

def affinity(r, s):
    return (r @ s + (1+math.sqrt(1-(r @ r)))(1+math.sqrt(1-(s @ s))))/((math.sqrt(1+LA.vector_norm(r))+math.sqrt(1-LA.vector_norm(r)))*(math.sqrt(1+LA.vector_norm(s))+math.sqrt(1-LA.vector_norm(s))))

def wooters_distance(r, s):
    return math.acos(math.sqrt(fidelity(r, s)))