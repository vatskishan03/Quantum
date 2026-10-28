#!/usr/bin/env python3
"""Second check, without stim: replays the faults in a transcription of decoder_setup.py from
github.com/sbravyi/BivariateBicycleCodes (commit fa77e33). Needs only python3 + numpy.

  python3 verify_transcription.py

Names follow decoder_setup.py (lin_order, nbs, cycle, sX, sZ, simulate_circuitZ / simulate_circuitX, syndrome
sparsification), with the [[288,12,18]] parameters that the file lists in a comment and num_cycles = 18.
Each fault is inserted into cycle_repeated the way decoder_setup.py builds its single-fault circuits
(head + [error] + tail, i.e. after the gate; before the gate for measurement errors).
For 1 noise-free cycle (Beverland et al.) and 2 (decoder_setup.py) it checks:
  (a) every sparsified syndrome bit is 0 (no detection event anywhere),
  (b) the final data error has zero syndrome for the checks that detect it,
  (c) the final data error is not a stabilizer, i.e. it is a nontrivial logical operator (GF(2) rank test, basis-free).
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ell, m = 12, 12                       # [[288,12,18]]
a1, a2, a3 = 3, 2, 7
b1, b2, b3 = 3, 1, 2
num_cycles = 18
n = 2 * m * ell
n2 = m * ell
I_ell = np.identity(ell, dtype=int); I_m = np.identity(m, dtype=int)
x = {i: np.kron(np.roll(I_ell, i, axis=1), I_m) for i in range(ell)}
y = {i: np.kron(I_ell, np.roll(I_m, i, axis=1)) for i in range(m)}
A = (x[a1] + y[a2] + y[a3]) % 2
B = (y[b1] + x[b2] + x[b3]) % 2
A1, A2, A3, B1, B2, B3 = x[a1], y[a2], y[a3], y[b1], x[b2], x[b3]
hx = np.hstack((A, B)); hz = np.hstack((B.T, A.T))

lin_order = {}; data_qubits = []; Xchecks = []; Zchecks = []
cnt = 0
for i in range(n2): Xchecks.append(('Xcheck', i)); lin_order[('Xcheck', i)] = cnt; cnt += 1
for i in range(n2): data_qubits.append(('data_left', i)); lin_order[('data_left', i)] = cnt; cnt += 1
for i in range(n2): data_qubits.append(('data_right', i)); lin_order[('data_right', i)] = cnt; cnt += 1
for i in range(n2): Zchecks.append(('Zcheck', i)); lin_order[('Zcheck', i)] = cnt; cnt += 1

nbs = {}
for i in range(n2):
    c = ('Xcheck', i)
    nbs[(c, 0)] = ('data_left', np.nonzero(A1[i, :])[0][0]); nbs[(c, 1)] = ('data_left', np.nonzero(A2[i, :])[0][0])
    nbs[(c, 2)] = ('data_left', np.nonzero(A3[i, :])[0][0]); nbs[(c, 3)] = ('data_right', np.nonzero(B1[i, :])[0][0])
    nbs[(c, 4)] = ('data_right', np.nonzero(B2[i, :])[0][0]); nbs[(c, 5)] = ('data_right', np.nonzero(B3[i, :])[0][0])
    c = ('Zcheck', i)
    nbs[(c, 0)] = ('data_left', np.nonzero(B1[:, i])[0][0]); nbs[(c, 1)] = ('data_left', np.nonzero(B2[:, i])[0][0])
    nbs[(c, 2)] = ('data_left', np.nonzero(B3[:, i])[0][0]); nbs[(c, 3)] = ('data_right', np.nonzero(A1[:, i])[0][0])
    nbs[(c, 4)] = ('data_right', np.nonzero(A2[:, i])[0][0]); nbs[(c, 5)] = ('data_right', np.nonzero(A3[:, i])[0][0])

sX = ['idle', 1, 4, 3, 5, 0, 2]
sZ = [3, 5, 0, 1, 2, 4, 'idle']
cycle = []
t = 0                                                           # round 0
for q in Xchecks: cycle.append(('PrepX', q))
cnoted = []
for target in Zchecks:
    control = nbs[(target, sZ[t])]; cnoted.append(control); cycle.append(('CNOT', control, target))
for q in data_qubits:
    if q not in cnoted: cycle.append(('IDLE', q))
for t in range(1, 6):                                           # rounds 1-5
    for control in Xchecks: cycle.append(('CNOT', control, nbs[(control, sX[t])]))
    for target in Zchecks: cycle.append(('CNOT', nbs[(target, sZ[t])], target))
t = 6                                                           # round 6
for q in Zchecks: cycle.append(('MeasZ', q))
cnoted = []
for control in Xchecks:
    target = nbs[(control, sX[t])]; cycle.append(('CNOT', control, target)); cnoted.append(target)
for q in data_qubits:
    if q not in cnoted: cycle.append(('IDLE', q))
for q in data_qubits: cycle.append(('IDLE', q))                # round 7
for q in Xchecks: cycle.append(('MeasX', q))
for q in Zchecks: cycle.append(('PrepZ', q))
cycle_repeated = num_cycles * cycle


def simulate_circuitZ(C):
    syndrome_history = []; syndrome_map = {}; state = np.zeros(2 * n, dtype=int); err_cnt = 0; syn_cnt = 0
    for gate in C:
        if gate[0] == 'CNOT':
            control = lin_order[gate[1]]; target = lin_order[gate[2]]
            state[control] = (state[target] + state[control]) % 2
            continue
        if gate[0] == 'PrepX':
            state[lin_order[gate[1]]] = 0; continue
        if gate[0] == 'MeasX':
            syndrome_history.append(state[lin_order[gate[1]]])
            syndrome_map.setdefault(gate[1], []).append(syn_cnt); syn_cnt += 1
            continue
        if gate[0] in ['Z', 'Y']:
            err_cnt += 1; q = lin_order[gate[1]]; state[q] = (state[q] + 1) % 2; continue
        if gate[0] in ['ZZ', 'YY', 'YZ', 'ZY']:
            err_cnt += 1
            for g in gate[1:]: q = lin_order[g]; state[q] = (state[q] + 1) % 2
            continue
    return np.array(syndrome_history, dtype=int), state, syndrome_map, err_cnt


def simulate_circuitX(C):
    syndrome_history = []; syndrome_map = {}; state = np.zeros(2 * n, dtype=int); err_cnt = 0; syn_cnt = 0
    for gate in C:
        if gate[0] == 'CNOT':
            control = lin_order[gate[1]]; target = lin_order[gate[2]]
            state[target] = (state[target] + state[control]) % 2
            continue
        if gate[0] == 'PrepZ':
            state[lin_order[gate[1]]] = 0; continue
        if gate[0] == 'MeasZ':
            syndrome_history.append(state[lin_order[gate[1]]])
            syndrome_map.setdefault(gate[1], []).append(syn_cnt); syn_cnt += 1
            continue
        if gate[0] in ['X', 'Y']:
            err_cnt += 1; q = lin_order[gate[1]]; state[q] = (state[q] + 1) % 2; continue
        if gate[0] in ['XX', 'YY', 'XY', 'YX']:
            err_cnt += 1
            for g in gate[1:]: q = lin_order[g]; state[q] = (state[q] + 1) % 2
            continue
    return np.array(syndrome_history, dtype=int), state, syndrome_map, err_cnt


def gf2_rank(M):
    basis = {}
    for row in M:
        r = int(''.join(str(int(v)) for v in row), 2)
        while r:
            h = r.bit_length() - 1
            if h in basis: r ^= basis[h]
            else: basis[h] = r; break
    return len(basis)


def check(sector, faults, noise_free):
    before, after = {}, {}
    for f in faults:
        b = f['bravyi']
        err = (b['error'][0],) + tuple((nm, i) for nm, i in b['error'][1:])
        gate = cycle_repeated[b['position_in_cycle_repeated']]
        assert [gate[0]] + [[q[0], int(q[1])] for q in gate[1:]] == b['gate'], f'gate mismatch at {b}'
        assert all(q in gate[1:] for q in err[1:]), 'error not on the qubits of its gate'
        # the fault must be an entry of decoder_setup.py's single-fault list for this sector
        P1, P2, meas, prep = ('Z', 'ZZ', 'MeasX', 'PrepX') if sector == 'Z' else ('X', 'XX', 'MeasZ', 'PrepZ')
        side = b['error_inserted']
        legal = ((gate[0] == meas and side == 'before gate' and err == (P1, gate[1])) or
                 (gate[0] in (prep, 'IDLE') and side == 'after gate' and err == (P1, gate[1])) or
                 (gate[0] == 'CNOT' and side == 'after gate' and err in ((P1, gate[1]), (P1, gate[2]), (P2, gate[1], gate[2]))))
        assert legal, f'not a single fault of decoder_setup.py: {b}'
        (before if b['error_inserted'] == 'before gate' else after).setdefault(b['position_in_cycle_repeated'], []).append(err)
    C = []
    for pos, g in enumerate(cycle_repeated):
        C += before.get(pos, []); C.append(g); C += after.get(pos, [])
    sim, checks = (simulate_circuitZ, Xchecks) if sector == 'Z' else (simulate_circuitX, Zchecks)
    syn, state, smap, err_cnt = sim(C + noise_free * cycle)
    s2 = syn.copy()                                             # syndrome sparsification as in decoder_setup.py
    for c in checks:
        pos = smap[c]
        for row in range(1, len(pos)): s2[pos[row]] = (syn[pos[row]] + syn[pos[row - 1]]) % 2
    e = np.array([state[lin_order[q]] for q in data_qubits])
    h_detect, h_same = (hx, hz) if sector == 'Z' else (hz, hx)
    a_ = not s2.any(); b_ = not ((h_detect @ e) % 2).any(); c_ = gf2_rank(np.vstack([h_same, e])) > gf2_rank(h_same)
    return err_cnt, a_, b_, c_, int(e.sum())


ok = True
for sector in 'ZX':
    w = json.load(open(os.path.join(HERE, f'witness_{sector}type_16faults.json')))
    for nf in (1, 2):
        err_cnt, a_, b_, c_, wt = check(sector, w['faults'], nf)
        good = err_cnt == 16 and a_ and b_ and c_
        ok &= good
        print(f'{sector}-type, {num_cycles} noisy + {nf} noise-free cycles: faults={err_cnt} no detection event={a_} '
              f'final data error weight={wt} zero code syndrome={b_} nontrivial logical={c_}  {"PASS" if good else "FAIL"}')
print('ALL CHECKS PASS' if ok else 'SOME CHECK FAILED')
