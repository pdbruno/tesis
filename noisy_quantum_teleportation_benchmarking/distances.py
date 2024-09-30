from numpy import linalg as LA
import numpy as np
from numpy.typing import NDArray
from typing import Callable


def trace_distance(r, s):
    return np.sqrt((r - s) @ (r - s)) / 2


def fidelity(r, s):
    return np.maximum(0, np.minimum(1, (1 + r @ s) / 2))


def affinity(r, s):
    return (
        r @ s
        + (1 + np.sqrt(1 - np.minimum(1, r @ r)))
        * (1 + np.sqrt(1 - np.minimum(1, s @ s)))
    ) / (
        (
            np.sqrt(1 + np.minimum(1, LA.vector_norm(r)))
            + np.sqrt(1 - np.minimum(1, LA.vector_norm(r)))
        )
        * (
            np.sqrt(1 + np.minimum(1, LA.vector_norm(s)))
            + np.sqrt(1 - np.minimum(1, LA.vector_norm(s)))
        )
    )


def wooters_distance(r, s):
    return np.acos(np.sqrt(fidelity(r, s)))


type DistanceFunction = Callable[
    [NDArray[np.floating], NDArray[np.floating]], NDArray[np.floating]
]
