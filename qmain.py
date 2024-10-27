from api_tokens import api_tokens
from noisy_quantum_teleportation_benchmarking.QPU_optimized.experiment_QPU_optimized import (
    ExperimentQPUOptimized,
)
from noisy_quantum_teleportation_benchmarking.experiment import Experiment
from noisy_quantum_teleportation_benchmarking.channels import AmplitudeDampingChannel
from noisy_quantum_teleportation_benchmarking.distances import fidelity
from noisy_quantum_teleportation_benchmarking.sampler import PauliSampler
import numpy as np
from qiskit.compiler import transpile
from qiskit_ibm_runtime import (
    QiskitRuntimeService,
    SamplerV2 as Sampler,
    IBMRuntimeError,
)

adc = AmplitudeDampingChannel()
channel_combinations = [(adc, adc)]
diagonal_exploration_space = [(x, x) for x in np.linspace(0, 1, 11)]
fighting_noise_with_noise_exploration_space = [(0.85, x) for x in np.linspace(0, 1, 11)]
distances = [np.vectorize(fidelity, signature="(3),(3)->()")]
input_sampler = PauliSampler()


class ServiceContainer:
    def __init__(self) -> None:
        self.iterador_tokens = iter(api_tokens)
        self.service: QiskitRuntimeService = None  # type: ignore
        self.log_in_con_proximo_token()

    def log_in_con_proximo_token(self):
        next_token = next(self.iterador_tokens)
        QiskitRuntimeService.save_account(
            channel="ibm_quantum",
            overwrite=True,
            token=api_tokens[next_token],
        )
        self.service = QiskitRuntimeService()
        print('el proximo token es de', next_token)


class QPUExperiments:
    def __init__(self, serviceContainer, backend_name) -> None:
        self.serviceContainer: ServiceContainer = serviceContainer
        self.backend_name = backend_name
        self.transpiler = lambda qc: transpile(
            qc,
            serviceContainer.service.backend(backend_name),
            optimization_level=3,
        )

    def try_run_or_change_account(self, fn):
        try:
            fn()
        except IBMRuntimeError as e:
            if "Job create exceeds open plan job usage limits" in e.args[0]:
                serviceContainer.log_in_con_proximo_token()
                fn()
            else:
                raise e

    def run_unoptimized_experiment(self, shots):

        experiment = Experiment(
            channel_combinations,
            diagonal_exploration_space,
            input_sampler,
            self.transpiler,
        )
        self.try_run_or_change_account(
            lambda: experiment.run_with_sampler(
                self.get_basic_sampler(),
                distances,
                shots,
                f"qpu-results\\{self.backend_name}_diagonal_unoptimized_{shots}.csv",
            )
        )

    def run_QPU_optimized_experiment(self, shots):

        experiment = ExperimentQPUOptimized(
            channel_combinations, diagonal_exploration_space, self.transpiler
        )
        self.try_run_or_change_account(
            lambda: experiment.run_with_sampler(
                self.get_basic_sampler(),
                distances,
                shots,
                f"qpu-results\\{self.backend_name}_diagonal_QPU_optimized_{shots}.csv",
            )
        )

    def run_QPU_optimized_EM_experiment(self, shots):
        experiment = ExperimentQPUOptimized(
            channel_combinations, diagonal_exploration_space, self.transpiler
        )
        self.try_run_or_change_account(
            lambda: experiment.run_with_sampler(
                self.get_EM_sampler(),
                distances,
                shots,
                f"qpu-results\\{self.backend_name}_diagonal_QPU_optimized_EM_{shots}.csv",
            )
        )

    def run_fighting_noise_with_noise_experiment(self, shots):
        experiment = ExperimentQPUOptimized(
            channel_combinations,
            fighting_noise_with_noise_exploration_space,
            self.transpiler,
        )
        self.try_run_or_change_account(
            lambda: experiment.run_with_sampler(
                self.get_EM_sampler(),
                distances,
                shots,
                f"qpu-results\\{self.backend_name}_fighting_noise_with_noise_{shots}.csv",
            )
        )

    def run_experiments(self):
        self.run_unoptimized_experiment(2000)
        self.run_QPU_optimized_experiment(2000)
        self.run_QPU_optimized_EM_experiment(2000)
        self.run_fighting_noise_with_noise_experiment(10000)

    def get_basic_sampler(self):
        return Sampler(mode=self.serviceContainer.service.backend(self.backend_name))

    def get_EM_sampler(self):
        sampler = Sampler(mode=self.serviceContainer.service.backend(self.backend_name))
        # Turn on gate twirling. Requires qiskit_ibm_runtime 0.23.0 or later.
        sampler.options.twirling.enable_gates = True  # type: ignore

        sampler.options.dynamical_decoupling.enable = True  # type: ignore
        sampler.options.dynamical_decoupling.sequence_type = "XpXm"  # type: ignore
        return sampler


serviceContainer = ServiceContainer()


for backend_name in ["ibm_sherbrook", "ibm_brisbane"]:#, "ibm_kyiv"
    QPUExperiments(serviceContainer, backend_name).run_experiments()

