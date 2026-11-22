# Cirac-Zoller CNOT with two trapped ions (Schmidt-Kaler et al., Nature 422, 408, 2003)
# Hilbert space: ion 1 x ion 2 x bus mode truncated at N phonons. |S> = 0, |D> = 1.
import numpy as np
from scipy.linalg import expm

N = 8
I2, IN = np.eye(2), np.eye(N)
sp = np.array([[0, 0], [1, 0]], complex)
sm = sp.conj().T
sx = sp + sm
sy = -1j * sp + 1j * sm
a = np.diag(np.sqrt(np.arange(1, N)), 1).astype(complex)
ad = a.conj().T
r2 = np.sqrt(2)
Xi = np.array([[0, 1], [1, 0]], complex)


def k3(A, B, C):
    return np.kron(np.kron(A, B), C)


def carrier_G(j, ph, eps=0.0):
    s = np.cos(ph) * sx + np.sin(ph) * sy
    t, o = (k3(s, I2, IN), k3(I2, s, IN)) if j == 1 else (k3(I2, s, IN), k3(s, I2, IN))
    return t + eps * o   # eps: part of the beam hitting the other ion


def blue_G(j, ph, eps=0.0):
    f = lambda A, B: np.exp(1j * ph) * k3(A, B, ad)
    t, o = (f(sp, I2), f(I2, sp)) if j == 1 else (f(I2, sp), f(sp, I2))
    G = t + eps * o
    return G + G.conj().T


# Table 1 of the paper: (carrier or blue sideband, ion, angle, phase)
SEQ = [("b", 1, np.pi, 0), ("c", 2, np.pi / 2, 0), ("b", 2, np.pi, 0), ("b", 2, np.pi / r2, np.pi / 2),
       ("b", 2, np.pi, 0), ("b", 2, np.pi / r2, np.pi / 2), ("c", 2, np.pi / 2, np.pi), ("b", 1, np.pi, np.pi)]


def pulse(kind, j, th, ph, eps=0.0, amp=1.0, dph=0.0):
    G = (carrier_G if kind == "c" else blue_G)(j, ph + dph, eps)
    return expm(-1j * amp * th / 2 * G)


def gate(eps=0.0, amp_sig=0.0, ph_sig=0.0, rng=None):
    # one run of the sequence, with a random amplitude and phase error on every pulse
    U = np.eye(4 * N, dtype=complex)
    for kind, j, th, ph in SEQ:
        amp = 1 + (rng.normal(0, amp_sig) if amp_sig else 0)
        dph = rng.normal(0, ph_sig) if ph_sig else 0
        U = pulse(kind, j, th, ph, eps, amp, dph) @ U
    return U


def vec(i1, i2, n=0):
    v = np.zeros(4 * N, complex)
    v[(i1 * 2 + i2) * N + n] = 1
    return v


def run(psi_in_fn, shots=200, nth=0.0, **err):
    # 4x4 density matrix of the two ions (bus traced out), averaged over shots
    rng = np.random.default_rng(1)
    rho = np.zeros((4, 4), complex)
    for _ in range(shots):
        n0 = 1 if rng.random() < nth else 0
        out = gate(rng=rng, **err) @ psi_in_fn(n0)
        M = out.reshape(4, N)
        rho += M @ M.conj().T
    return rho / shots


def truth_table(**kw):
    T = np.zeros((4, 4))
    for i, (x, y) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
        T[i] = np.real(np.diag(run(lambda n0: vec(x, y, n0), **kw)))
    return T


def diag_mean(T):
    return np.mean([T[0, 0], T[1, 1], T[2, 3], T[3, 2]])


def bell_input(n0):
    return (vec(0, 0, n0) + vec(1, 0, n0)) / np.sqrt(2)


def bell_fidelity(rho):
    # target (|SS> + e^{i alpha}|DD>)/sqrt(2), best alpha
    best = 0
    for al in np.linspace(0, 2 * np.pi, 361):
        tgt = np.zeros(4, complex)
        tgt[0], tgt[3] = 1 / r2, np.exp(1j * al) / r2
        best = max(best, np.real(tgt.conj() @ rho @ tgt))
    return best


def parity_curve(rho, phis):
    # pi/2 pulse on both ions with phase phi, parity = P(SS) + P(DD) - P(SD) - P(DS)
    out = []
    for phi in phis:
        R = expm(-1j * np.pi / 4 * (np.cos(phi) * sx + np.sin(phi) * sy))
        R2 = np.kron(R, R)
        p = np.real(np.diag(R2 @ rho @ R2.conj().T))
        out.append(p[0] + p[3] - p[1] - p[2])
    return np.array(out)


# --- time-correlated laser noise ---
OMEGA_B = np.pi / 95.0   # rad/us, blue sideband pi pulse = 95 us (paper)
OMEGA_C = 1.2566         # rad/us, carrier pi/2 ~ 1.25 us


def duration(kind, th):
    return th / (OMEGA_B if kind == "b" else OMEGA_C)


GATE_T = sum(duration(k, th) for k, j, th, ph in SEQ)


def laser_phase(rng, T, sig_d, tau, dt=1.0):
    # laser phase = integral of an Ornstein-Uhlenbeck frequency noise (std sig_d rad/us, correlation time tau us)
    n = int(T / dt) + 2
    if tau <= 0:
        d = rng.normal(0, sig_d, n)
    else:
        d = np.zeros(n)
        c = np.exp(-dt / tau)
        d[0] = rng.normal(0, sig_d)
        for i in range(1, n):
            d[i] = c * d[i - 1] + rng.normal(0, sig_d * np.sqrt(1 - c * c))
    return np.cumsum(d) * dt


def sigma_for_phase(Phi, tau, T=GATE_T):
    # frequency noise std that gives an accumulated phase std of Phi over the gate
    if tau <= 0:
        return Phi / np.sqrt(T)
    return Phi / np.sqrt(2 * tau ** 2 * (T / tau - 1 + np.exp(-T / tau)))


def gate_corr(rng, sig_d, tau, eps=0.0):
    # each pulse gets the laser phase at its midpoint
    durs = [duration(k, th) for k, j, th, ph in SEQ]
    phiL = laser_phase(rng, sum(durs), sig_d, tau)
    t, U = 0.0, np.eye(4 * N, dtype=complex)
    for (kind, j, th, ph), d in zip(SEQ, durs):
        U = pulse(kind, j, th, ph, eps, 1.0, phiL[int(t + d / 2)]) @ U
        t += d
    return U


# --- deterministic versions, bus temperature ---
def gate_amp(e):
    U = np.eye(4 * N, dtype=complex)
    for kind, j, th, ph in SEQ:
        U = pulse(kind, j, th, ph, 0.0, 1 + e, 0.0) @ U
    return U


def bell_rho(U, n0=0):
    M = (U @ bell_input(n0)).reshape(4, N)
    return M @ M.conj().T


def thermal_weights(nbar):
    p = nbar ** np.arange(N) / (1 + nbar) ** (np.arange(N) + 1) if nbar > 0 else np.eye(N)[0]
    return p / p.sum()


def eff_unitary(U):
    # 4x4 block of the ions when the bus starts and ends in |n=0>
    return U.reshape(2, 2, N, 2, 2, N)[:, :, 0, :, :, 0].reshape(4, 4)


# --- Molmer-Sorensen gate (bichromatic, Lamb-Dicke limit) ---
def ms_gate(g, delta=2 * np.pi / 300.0, steps=300):
    # H(t) = g Jx (a e^{i delta t} + h.c.) over T = 2 pi / delta
    J = k3(Xi, I2, IN) + k3(I2, Xi, IN)
    A, Ad = k3(I2, I2, a), k3(I2, I2, ad)
    dt = 2 * np.pi / delta / steps
    U = np.eye(4 * N, dtype=complex)
    for k in range(steps):
        t = (k + 0.5) * dt
        U = expm(-1j * g * J @ (A * np.exp(1j * delta * t) + Ad * np.exp(-1j * delta * t)) * dt) @ U
    return U


def ms_ideal():
    return expm(-1j * np.pi / 4 * np.kron(Xi, Xi))


def ms_bell_rho(U, n0=0):
    M = (U @ vec(0, 0, n0)).reshape(4, N)
    return M @ M.conj().T
