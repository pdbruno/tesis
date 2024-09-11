import numpy as np

bell_states = ((0, 0), (0, 1), (1, 0), (1, 1))


def get_counts_for_bell_state(
    bell_state: tuple[int, int], measurements_x, measurements_y, measurements_z
):
    return (
        measurements_x.get(bell_state + (0,), 0)
        + measurements_x.get(bell_state + (1,), 0)
        + measurements_y.get(bell_state + (0,), 0)
        + measurements_y.get(bell_state + (1,), 0)
        + measurements_z.get(bell_state + (0,), 0)
        + measurements_z.get(bell_state + (1,), 0)
    )


def get_coordinate(base, bell_state):
    zeroes = base.get(bell_state + (0,), 0)
    ones = base.get(bell_state + (1,), 0)
    return (zeroes - ones) / (zeroes + ones) if zeroes + ones != 0 else 0

def get_probabilities_and_states(measurements_xyz, shots):
    measurements_x, measurements_y, measurements_z = measurements_xyz
    probabilities = np.array(
        [
            get_counts_for_bell_state(
                bell_state, measurements_x, measurements_y, measurements_z
            )
            / (shots * 3)
            for bell_state in bell_states
        ]
    )
    states = np.array(
        [
            [
                get_coordinate(measurements_x, bell_state),
                get_coordinate(measurements_y, bell_state),
                get_coordinate(measurements_z, bell_state),
            ]
            for bell_state in bell_states
        ]
    )
    return probabilities, states

def get_counts(result_for_basis_measurement):
    bitstrings, counts = np.unique(
        result_for_basis_measurement.T,
        axis=0,
        return_counts=True,
    )
    return {
        tuple(bitstring.data): count.item()
        for bitstring, count in zip(bitstrings, counts)
    }

def quantum_state_tomography(result_for_input_state, shots):
    return get_probabilities_and_states(
        [get_counts(result_for_basis_measurement) for result_for_basis_measurement in result_for_input_state],
        shots,
    )