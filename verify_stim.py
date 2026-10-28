#!/usr/bin/env python3
"""Check the 16-fault logical errors of IBM's two-gross memory circuit with stim (needs: python3, numpy, stim).

  python3 verify_stim.py

1. Model check: builds the noisy memory circuits (gen_bb_memory.py) and checks that they have the sizes Beverland et al.
   (arXiv:2511.15177, Table 1) list for the gross and two-gross circuits under Z-type noise (detectors M, observables K,
   total fault weight N in units of p/15, distinct fault signatures N~). Matching sizes is a necessary condition only; the
   script also prints an exact fingerprint of the two-gross Z-type system (sha256 of its sorted fault signatures with
   multiplicities, see README) that can be compared with one computed from decoder_setup.py's own fault list.
2. Witnesses: inserts each listed fault as a probability-1 Pauli error right after its CNOT into the noise-free circuit,
   after checking that it is a single CNOT fault of a noisy cycle (the noisy circuit applies DEPOLARIZE2 to exactly that
   CNOT pair, and the error acts only on that pair). Stim then samples the faulty circuit: a logical error with no
   detection event at all means zero detectors fire and at least one logical observable flips.
3. The same check for a term-order variant of the schedule (roles of A2 = y^2 and A3 = y^7 exchanged).
"""
import json, os, sys
import numpy as np
import stim

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gen_bb_memory import build

BEVERLAND_TABLE1 = {  # QEC system: (N~, N, M, K) as printed in arXiv:2511.15177 Table 1
    'gross': ('BB(12)-circuit', 8785, 184896, 936, 12, 12),
    'two-gross': ('BB(18)-circuit', 26209, 554688, 2736, 12, 18),
}
UNITS = {  # weight, in units of p/15, of the Z-type (basis X) or X-type (basis Z) part of each noise channel of strength p
    'X': {'Z_ERROR': 15, 'MX': 15, 'DEPOLARIZE1': 10, 'DEPOLARIZE2': 12},
    'Z': {'X_ERROR': 15, 'M': 15, 'DEPOLARIZE1': 10, 'DEPOLARIZE2': 12},
}


def system_sizes(circ, basis):
    dem = circ.detector_error_model(decompose_errors=False)      # also checks that detectors/observables are deterministic
    sigs = set()
    for ins in dem.flattened():
        if ins.type != 'error': continue
        key = tuple(sorted((t.is_logical_observable_id(), t.val) for t in ins.targets_copy()))
        if any(not o for o, _ in key): sigs.add(key)            # signatures that fire at least one detector
    N = 0
    for ins in circ:
        u = UNITS[basis].get(ins.name)
        args = ins.gate_args_copy()
        if u and args and args[0] > 0:
            n = len(ins.targets_copy())
            N += u * (n // 2 if ins.name == 'DEPOLARIZE2' else n)
    return circ.num_detectors, circ.num_observables, N, len(sigs)


def fingerprint(code, R):
    """sha256 of the Z-type system: one line per distinct nonzero detector signature, 'd1,d2,...:u' (detectors in the
    order cycle-major, X checks in decoder_setup.py order; u = fault weight in units of p/15), sorted, then 'zero:u0'."""
    import hashlib
    p = 1e-9
    circ = stim.Circuit(build(code, 'X', R, p)[0])
    items, tot = [], 0
    for ins in circ.detector_error_model(decompose_errors=False).flattened():
        if ins.type != 'error': continue
        ds = sorted(t.val for t in ins.targets_copy() if t.is_relative_detector_id())
        if not ds: continue
        u = round(ins.args_copy()[0] / (p / 15)); tot += u
        items.append(','.join(map(str, ds)) + ':' + str(u))
    N = system_sizes(stim.Circuit(build(code, 'X', R, 0.001)[0]), 'X')[2]
    items.sort(); items.append(f'zero:{N - tot}')
    return hashlib.sha256('\n'.join(items).encode()).hexdigest(), len(items) - 1, N - tot


def inject(noiseless, noisy, faults, R):
    """noise-free circuit + one probability-1 Pauli error after the CNOT of each fault; every fault is checked to be a
    single fault of the noise model (on one CNOT pair that the noisy circuit hits with DEPOLARIZE2, in a noisy cycle)."""
    by_layer = {}
    for f in faults: by_layer.setdefault(f['cx_layer'], []).append(f)
    noisy_layers = []                       # (CX targets, following instruction) in the noisy circuit
    insts = list(noisy)
    for i, ins in enumerate(insts):
        if ins.name == 'CX': noisy_layers.append((ins, insts[i + 1] if i + 1 < len(insts) else None))
    out = stim.Circuit(); k = -1; placed = 0
    for ins in noiseless:
        out.append(ins)
        if ins.name != 'CX': continue
        k += 1
        if k not in by_layer: continue
        t = [x.value for x in ins.targets_copy()]
        pairs = set(zip(t[0::2], t[1::2]))
        ncx, nxt = noisy_layers[k]
        assert [x.value for x in ncx.targets_copy()] == t, 'noise-free and noisy circuits differ at this CNOT layer'
        assert nxt is not None and nxt.name == 'DEPOLARIZE2' and [x.value for x in nxt.targets_copy()] == t, \
            f'CNOT layer {k} is not followed by DEPOLARIZE2 on the same pairs in the noisy circuit'
        for f in by_layer[k]:
            c, tg = f['cnot']; qs = f['qubits']; P = f['pauli'][0]
            assert (c, tg) in pairs, f'fault {f} is not on a CNOT of layer {k}'
            assert set(qs) <= {c, tg} and len(qs) in (1, 2) and f['cycle'] < R and k == 7 * f['cycle'] + f['round']
            if len(qs) == 1: out.append(f'{P}_ERROR', qs, 1.0)
            else: out.append('CORRELATED_ERROR', [stim.target_z(q) if P == 'Z' else stim.target_x(q) for q in qs], 1.0)
            placed += 1
    assert placed == len(faults), 'some faults were not placed'
    return out


def check_witness(code, sector, basis, fname, tag, R=18):
    w = json.load(open(os.path.join(HERE, fname)))
    noiseless = stim.Circuit(build(code, basis, R, 0.001, noiseless=True)[0])
    noisy = stim.Circuit(build(code, basis, R, 0.001)[0])
    assert noisy.without_noise() == noiseless, 'noise-free circuit is not the noisy circuit with noise removed'
    faulty = inject(noiseless, noisy, w['faults'], R)
    dets, obs = faulty.compile_detector_sampler().sample(shots=1, separate_observables=True)
    d0, o0 = noiseless.compile_detector_sampler().sample(shots=1, separate_observables=True)
    nf = len(w['faults']); fired = int(dets.sum()); flips = [k for k in range(obs.shape[1]) if obs[0, k]]
    # minimality: detection events are linear in the faults, so check all 2^16 subsets from the 16 single-fault patterns
    single = []
    for f in w['faults']:
        dd, _ = inject(noiseless, noisy, [f], R).compile_detector_sampler().sample(shots=1, separate_observables=True)
        single.append(int(''.join('1' if b else '0' for b in dd[0]), 2))
    acc = [0] * (1 << nf)
    for s in range(1, 1 << nf):
        low = s & -s; acc[s] = acc[s ^ low] ^ single[low.bit_length() - 1]
    silent = [s for s in range(1, 1 << nf) if acc[s] == 0]
    minimal = silent == [(1 << nf) - 1]
    good = nf == 16 and fired == 0 and len(flips) > 0 and not d0.any() and not o0.any() and minimal
    print(f'{tag} {sector}-type faults ({basis}-basis memory, {R} noisy cycles + 1 noise-free): {nf} single CNOT faults, '
          f'detectors fired: {fired} of {faulty.num_detectors}, logical observables flipped: {flips}, '
          f'no proper subset undetected: {minimal}  {"PASS" if good else "FAIL"}')
    return good


def main():
    ok = True
    print('1. Consistency with Beverland et al. Table 1 (Z-type noise, standard circuit noise, R = d noisy cycles + 1 noise-free cycle;'
          ' a necessary condition, see README)')
    for code, (label, Nt, N_, M_, K_, R) in BEVERLAND_TABLE1.items():
        circ = stim.Circuit(build(code, 'X', R, 0.001)[0])
        M, K, N, sig = system_sizes(circ, 'X')
        good = (M, K, N, sig + 1) == (M_, K_, N_, Nt)
        ok &= good
        print(f'   {label:15s} detectors M={M} (paper {M_})  observables K={K} (paper {K_})  fault weight N={N} (paper {N_})  '
              f'distinct nonzero signatures {sig} (+1 all-zero column = {sig + 1}; paper N~={Nt})  {"MATCH" if good else "MISMATCH"}')
    fp, nsig, zero = fingerprint('two-gross', 18)
    print(f'   two-gross Z-type system fingerprint: sha256 {fp} ({nsig} signatures + zero class of weight {zero})')
    for sector, basis in (('Z', 'X'), ('X', 'Z')):
        ok &= check_witness('two-gross', sector, basis, f'witness_{sector}type_16faults.json', '2.')
    print('3. Term-order variant: the same code with the roles of A2 = y^2 and A3 = y^7 exchanged in the schedule')
    for sector, basis in (('Z', 'X'), ('X', 'Z')):
        ok &= check_witness('two-gross-A23swap', sector, basis, f'witness_A23swap_{sector}type_16faults.json', '  ')
    print('ALL CHECKS PASS' if ok else 'SOME CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
