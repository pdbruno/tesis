from abc import ABC, abstractmethod
from math import pi
import numpy as np
from numpy.typing import NDArray
from scipy.stats import rv_continuous


class BaseInputSampler(ABC):
    @abstractmethod
    def get_samples(self) -> NDArray[np.floating]:
        pass

    def __init__(self, length) -> None:
        super().__init__()
        self.length = length


class PauliSampler(BaseInputSampler):
    def __init__(self) -> None:
        super().__init__(6)

    eigenvector_x_mas = [pi / 2, 0, pi]  # [1., 0., 0.]
    eigenvector_x_menos = [-pi / 2, 0, pi]  # [-1., 0., 0.]
    eigenvector_y_mas = [-pi / 2, -pi / 2, pi / 2]  # [0., 1., 0.]
    eigenvector_y_menos = [pi / 2, -pi / 2, pi / 2]  # [0., -1., 0.]
    eigenvector_z_mas = [0.0, 0.0, 0.0]  # [0., 0., 1.]
    eigenvector_z_menos = [pi, 0, 0]  # [0., 0., -1.]

    all = np.array(
        [
            eigenvector_x_mas,
            eigenvector_x_menos,
            eigenvector_y_mas,
            eigenvector_y_menos,
            eigenvector_z_mas,
            eigenvector_z_menos,
        ]
    )

    def get_samples(self) -> NDArray[np.floating]:
        return self.all


# class sin_prob_dist(rv_continuous):
#     def _pdf(self, theta):
#         # The 0.5 is so that the distribution is normalized
#         return 0.5 * np.sin(theta)


class HaarMeasureSampler(BaseInputSampler):
    def __init__(self, length) -> None:
        super().__init__(length)
        # Samples of theta should be drawn from between 0 and pi
        # self.sin_sampler = sin_prob_dist(a=0, b=np.pi)

    def get_samples(self) -> NDArray[np.floating]:
        # thetas: NDArray[np.floating] = self.sin_sampler.rvs(size=self.length) # type: ignore
        thetas = np.arccos(-np.random.uniform(-1, 1, (self.length, 1)))
        phis_lambdas = np.random.uniform(0, 2 * pi, (self.length, 2))
        return np.hstack([thetas, phis_lambdas])
