#!/usr/bin/env python3
"""Stim circuits for IBM's bivariate-bicycle (BB) quantum memory experiment.

Follows decoder_setup.py of github.com/sbravyi/BivariateBicycleCodes (commit fa77e33) gate for gate: the depth-8
syndrome cycle with CNOT orders sX = [idle,1,4,3,5,0,2] and sZ = [3,5,0,1,2,4,idle], qubit numbering = its lin_order.
Uses the noise convention of Beverland et al., arXiv:2511.15177, Sec. 2.4: standard circuit noise with strength p on
R cycles, followed by ONE noise-free cycle; the logical action is read from the final data state.

Fault locations (as in decoder_setup.py): state preparation of check qubits (error after), check measurements (flip),
data qubits that are idle in a round (error after), every CNOT (two-qubit Pauli error after). No noise on data
preparation or final data measurement, none in the last cycle. The Z-check reset at the end of the last noisy cycle
prepares the noise-free cycle and is noise-free too; otherwise the X-type system gets artificial 4-fault logical errors
at the final boundary (a late data X error hidden from the noise-free cycle by three faulty resets). This choice does not
change the Z-type system.

  python3 gen_bb_memory.py CODE BASIS R OUT.stim [--p P] [--noiseless]
    CODE  = two-gross ([[288,12,18]]: ell=m=12, A=x^3+y^2+y^7, B=y^3+x+x^2) or gross ([[144,12,12]]: ell=12, m=6, A=x^3+y+y^2, B=y^3+x+x^2)
            or two-gross-A23swap (two-gross code, schedule with the roles of A2 = y^2 and A3 = y^7 exchanged)
    BASIS = X (data in |+>, X-check detectors: the system of Z-type faults, Beverland's H_Z) or Z (data in |0>, Z-check detectors)
    R     = number of noisy syndrome cycles (Beverland: R = code distance); one noise-free cycle is appended.

Only numpy and stim are needed. Detector coordinates are (a, b, cycle, type) with check index a*m + b, type 0 = X, 3 = Z.
"""
import sys
import numpy as np

CODES = {
    'two-gross': (12, 12, 3, 2, 7, 3, 1, 2),   # ell, m, a1, a2, a3, b1, b2, b3 (decoder_setup.py, "[[288,12,18]]")
    'gross': (12, 6, 3, 1, 2, 3, 1, 2),        # decoder_setup.py, "[[144,12,12]]"
    'two-gross-A23swap': (12, 12, 3, 7, 2, 3, 1, 2),   # same code, roles of A2 = y^2 and A3 = y^7 exchanged in the schedule
}


def logical_basis(h_perp, h_same):
    """Rows spanning ker(h_perp) modulo rowspace(h_same), over GF(2): logical operators of the measured type."""
    nq = h_perp.shape[1]
    to_int = lambda row: sum(int(v) << j for j, v in enumerate(row))
    pivots, R = [], []
    for r in map(to_int, h_perp):
        for pc, pr in zip(pivots, R):
            if (r >> pc) & 1: r ^= pr
        if r:
            pc = r.bit_length() - 1
            for j in range(len(R)):
                if (R[j] >> pc) & 1: R[j] ^= r
            pivots.append(pc); R.append(r)
    pset = set(pivots)
    kernel = []
    for f in range(nq):
        if f in pset: continue
        v = 1 << f
        for pc, pr in zip(pivots, R):
            if (pr >> f) & 1: v |= 1 << pc
        kernel.append(v)
    basis = {}
    for r in map(to_int, h_same):
        while r:
            h = r.bit_length() - 1
            if h in basis: r ^= basis[h]
            else: basis[h] = r; break
    out = []
    for v in kernel:
        w = v
        while w:
            h = w.bit_length() - 1
            if h in basis: w ^= basis[h]
            else: basis[h] = w; out.append(v); break
    return np.array([[(v >> j) & 1 for j in range(nq)] for v in out], dtype=np.uint8)


def build(code, basis, R, p, noiseless=False):
    ell, m, a1, a2, a3, b1, b2, b3 = CODES[code]
    n2 = ell * m
    Iell = np.identity(ell, dtype=int); Im = np.identity(m, dtype=int)
    x = {i: np.kron(np.roll(Iell, i, axis=1), Im) for i in range(ell)}
    y = {i: np.kron(Iell, np.roll(Im, i, axis=1)) for i in range(m)}
    A1, A2, A3, B1, B2, B3 = x[a1], y[a2], y[a3], y[b1], x[b2], x[b3]
    A = (A1 + A2 + A3) % 2; B = (B1 + B2 + B3) % 2
    hx = np.hstack((A, B)); hz = np.hstack((B.T, A.T))
    Xc = lambda i: i                 # lin_order: X checks, left data, right data, Z checks
    dL = lambda i: n2 + i
    dR = lambda i: 2 * n2 + i
    Zc = lambda i: 3 * n2 + i
    nz = lambda M_, i: int(np.nonzero(M_[i])[0][0])
    nb = {}
    for i in range(n2):                                   # neighbours, numbered 0..5 as in decoder_setup.py
        nb['X', i, 0], nb['X', i, 1], nb['X', i, 2] = dL(nz(A1, i)), dL(nz(A2, i)), dL(nz(A3, i))
        nb['X', i, 3], nb['X', i, 4], nb['X', i, 5] = dR(nz(B1, i)), dR(nz(B2, i)), dR(nz(B3, i))
        nb['Z', i, 0], nb['Z', i, 1], nb['Z', i, 2] = dL(nz(B1.T, i)), dL(nz(B2.T, i)), dL(nz(B3.T, i))
        nb['Z', i, 3], nb['Z', i, 4], nb['Z', i, 5] = dR(nz(A1.T, i)), dR(nz(A2.T, i)), dR(nz(A3.T, i))
    sX = ['idle', 1, 4, 3, 5, 0, 2]
    sZ = [3, 5, 0, 1, 2, 4, 'idle']
    data = [dL(i) for i in range(n2)] + [dR(i) for i in range(n2)]
    xs = ' '.join(str(Xc(i)) for i in range(n2)); zs = ' '.join(str(Zc(i)) for i in range(n2))
    L = []
    for i in range(n2):
        a, b = divmod(i, m)
        for q in (Xc(i), dL(i), dR(i), Zc(i)):
            L.append(f'QUBIT_COORDS({a}, {b}, {(q // n2)}) {q}')
    L.append(('RX ' if basis == 'X' else 'R ') + ' '.join(map(str, data)))
    L.append('R ' + zs)
    L.append('TICK')
    nmeas = 0
    rec = {}
    for c in range(R + 1):
        P = 0.0 if (noiseless or c == R) else p
        def cx(pairs):
            L.append('CX ' + ' '.join(f'{u} {v}' for u, v in pairs))
            if P: L.append(f'DEPOLARIZE2({P}) ' + ' '.join(f'{u} {v}' for u, v in pairs))
        def idle(qs):
            if P and qs: L.append(f'DEPOLARIZE1({P}) ' + ' '.join(map(str, qs)))
        # round 0: prepare X checks; CNOTs of the Z checks; the other data qubits idle
        L.append('RX ' + xs)
        if P: L.append(f'Z_ERROR({P}) ' + xs)
        pairs = [(nb['Z', i, sZ[0]], Zc(i)) for i in range(n2)]
        cx(pairs); used = {u for u, _ in pairs}; idle([q for q in data if q not in used])
        L.append('TICK')
        # rounds 1-5: CNOTs of the X checks (check = control) and of the Z checks (check = target)
        for t in range(1, 6):
            cx([(Xc(i), nb['X', i, sX[t]]) for i in range(n2)] + [(nb['Z', i, sZ[t]], Zc(i)) for i in range(n2)])
            L.append('TICK')
        # round 6: measure Z checks; last CNOTs of the X checks; the other data qubits idle
        L.append((f'M({P}) ' if P else 'M ') + zs)
        for i in range(n2): rec['Z', c, i] = nmeas + i
        nmeas += n2
        pairs = [(Xc(i), nb['X', i, sX[6]]) for i in range(n2)]
        cx(pairs); used = {v for _, v in pairs}; idle([q for q in data if q not in used])
        L.append('TICK')
        # round 7: all data idle; measure X checks; prepare Z checks
        idle(data)
        L.append((f'MX({P}) ' if P else 'MX ') + xs)
        for i in range(n2): rec['X', c, i] = nmeas + i
        nmeas += n2
        L.append('R ' + zs)
        if P and c < R - 1: L.append(f'X_ERROR({P}) ' + zs)   # the reset that feeds the noise-free cycle belongs to it
        L.append('TICK')
        for i in range(n2):
            a, b = divmod(i, m)
            cur = rec[basis, c, i] - nmeas
            ty = 0 if basis == 'X' else 3
            if c == 0: L.append(f'DETECTOR({a}, {b}, {c}, {ty}) rec[{cur}]')
            else: L.append(f'DETECTOR({a}, {b}, {c}, {ty}) rec[{cur}] rec[{rec[basis, c - 1, i] - nmeas}]')
    L.append(('MX ' if basis == 'X' else 'M ') + ' '.join(map(str, data)))
    drec = {q: nmeas + j for j, q in enumerate(data)}
    nmeas += len(data)
    lg = logical_basis(hz, hx) if basis == 'X' else logical_basis(hx, hz)
    for k, row in enumerate(lg):
        L.append(f'OBSERVABLE_INCLUDE({k}) ' + ' '.join(f'rec[{drec[data[j]] - nmeas}]' for j in np.nonzero(row)[0]))
    return '\n'.join(L) + '\n', len(lg)


if __name__ == '__main__':
    a = sys.argv[1:]
    code, basis, R, out = a[0], a[1], int(a[2]), a[3]
    p = float(a[a.index('--p') + 1]) if '--p' in a else 0.001
    txt, k = build(code, basis, R, p, noiseless='--noiseless' in a)
    open(out, 'w').write(txt)
    print(f'wrote {out}: {code} basis={basis} R={R} k={k} {"noiseless" if "--noiseless" in a else f"p={p}"}')
