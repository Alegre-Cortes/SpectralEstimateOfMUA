import numpy as np


def compute_moving_average_of_wf(wf, window_size, step_size):
    """
    Computes the moving average of a signal.

    Args:
        wf (dict): A dictionary containing the signal information.
            - wf["value"] (np.ndarray): The signal values.
            - wf["time"] (np.ndarray, optional): The corresponding time points (default: None).
        window_size (int): The size of the moving average window.
        step_size (int): The step size for moving the window.

    Returns:
        tuple:
            - t (np.ndarray): Time vector for the moving average.
            - ma (np.ndarray): Moving average of the signal.
    """

    num_samples = int(np.floor((len(wf["value"]) - window_size) / step_size)) + 1

    if "time" not in wf:
        dt = 1  # Default time step if not provided
    else:
        dt = wf["time"][1] - wf["time"][0]

    t = np.arange(num_samples) * dt * step_size + np.mean(wf["time"][:window_size])

    data = np.zeros((window_size, num_samples))
    for i in range(window_size):
        data[i, :] = wf["value"][i::step_size]

    ma = np.mean(data, axis=0)

    return t, ma
