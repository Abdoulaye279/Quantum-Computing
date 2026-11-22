import numpy as np
from cirac_zoller import (truth_table, run, bell_input, bell_fidelity, gate_amp, eff_unitary,
                          ms_gate, ms_ideal, sigma_for_phase, laser_phase, GATE_T)


def test_ideal_gate_is_a_cnot_in_populations():
    T = truth_table(shots=1)
    assert np.allclose(T, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], atol=1e-9)
    assert bell_fidelity(run(bell_input, shots=1)) > 0.999


def test_full_operator_is_minus_cy():
    # CY with ion 1 as control, ion 2 as target, global sign -1
    CY = np.eye(4, dtype=complex)
    CY[2:, 2:] = [[0, -1j], [1j, 0]]
    assert np.allclose(eff_unitary(gate_amp(0.0)), -CY, atol=1e-9)


def test_ms_gate_matches_ideal_xx():
    Ue = eff_unitary(ms_gate(2 * np.pi / 300.0 / 4))
    assert abs(np.trace(ms_ideal().conj().T @ Ue)) ** 2 / 16 > 0.999


def test_accumulated_phase_has_the_requested_std():
    for tau in [1, 100, 6000]:
        sd = sigma_for_phase(0.5, tau)
        end = [laser_phase(np.random.default_rng(s), GATE_T, sd, tau)[int(GATE_T)] for s in range(400)]
        assert abs(np.std(end) - 0.5) < 0.07
