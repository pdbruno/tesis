from numpy import linalg as LA
import numpy as np
from numpy.typing import NDArray
from typing import Callable


def vector_norm(x):
    return np.minimum(np.sqrt(x.dot(x)), 1)


def trace_distance(r, s):
    return vector_norm(r - s) / 2


def fidelity(r, s):
    return np.clip((1 + r.dot(s)) / 2, 0, 1)


def affinity(r, s):
    return (r.dot(s) + 1 + np.sqrt(1 - np.minimum(1, s.dot(s)))) / (
        np.sqrt(2) * (np.sqrt(1 + vector_norm(s)) + np.sqrt(1 - vector_norm(s)))
    )


def wooters_distance(r, s):
    return np.acos(np.sqrt(fidelity(r, s)))


type DistanceFunction = Callable[
    [NDArray[np.floating], NDArray[np.floating]], NDArray[np.floating]
]
