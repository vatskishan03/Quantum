# Verification details

## Result and circuit

Each witness is a list of 16 Pauli faults, inserted immediately after specified CNOT gates. The combined faults produce no detection events and act nontrivially on the encoded information. They therefore certify `d_circ <= 16` for the corresponding memory experiment.

The original schedule follows [`decoder_setup.py`](https://github.com/sbravyi/BivariateBicycleCodes/blob/fa77e3333d3ec44c79d8f914dd24c040d1da471b/decoder_setup.py) at commit `fa77e3333d3ec44c79d8f914dd24c040d1da471b`, with:

```text
ell = m = 12
(A1, A2, A3) = (x^3, y^2, y^7)
(B1, B2, B3) = (y^3, x, x^2)
sX = [idle, 1, 4, 3, 5, 0, 2]
sZ = [3, 5, 0, 1, 2, 4, idle]
```

There are 288 data qubits and 288 check ancillas. A syndrome cycle has eight rounds, numbered 0–7, including seven CNOT layers. The Stim checks use 18 noisy cycles followed by one noise-free cycle. Z-type faults are checked in X-basis memory; X-type faults in Z-basis memory.

The `A23swap` variant exchanges `A2 = y^2` and `A3 = y^7` in the schedule. It preserves the code while changing the gate order. Two separate witnesses certify the upper bound for this variant. Its relation to the schedule depicted in Figure 4b of [*Tour de gross*](https://arxiv.org/abs/2506.03094) is an interpretation of that figure, not an automatic comparison with the paper's supplied circuit files.

## Noise and final boundary

The generator places preparation errors after check preparation, measurement errors at check measurement, depolarizing noise on idle data qubits, and two-qubit depolarizing noise after CNOTs. Initial data preparation, final data measurement, and the closing cycle are noise-free.

The Z-check reset at the end of the last noisy cycle prepares the noise-free cycle and is treated as noise-free. This convention matters for the X-type fault system; it does not change the Z-type system used in the published size comparison. Claims for a different boundary convention require a separate model check.

For witness verification, stochastic noise is removed and the listed Pauli faults are applied with probability one. The sampler checks a fixed fault pattern, not a Monte Carlo estimate of its frequency. The noisy circuit is also built to check that every listed fault is an allowed fault at the specified CNOT.

## Two verification methods

`verify_stim.py` uses Stim to check all four witnesses. It checks the fault-free baseline, zero detection events, and a nonzero logical-observable change. It also enumerates all `2^16` subsets of each witness using linearity of the detector response. Only the full nonempty set is undetected. This is minimality within that particular set of faults, not a global distance lower bound.

`verify_transcription.py` replays the two original-schedule witnesses using NumPy and a transcription of the upstream circuit and Pauli propagation rules. It checks both one and two noise-free closing cycles. A binary rank test verifies that the final data error lies outside the stabilizer span, without using the Stim generator's logical-observable basis. In each case, the final data error has weight 24.

These are separate implementations, but share the underlying circuit specification and supplied witnesses. The transcription checker does not check the `A23swap` variant.

Use `python verify.py` for automated verification. It checks file integrity, requires both checkers to end with `ALL CHECKS PASS`, checks their process exit statuses, and compares the printed fingerprint with the value below. It runs the checkers with Python assertions enabled. The original transcription script alone prints its result but does not set a failing exit status for every failed check.

## Published system sizes

The Stim checker compares the Z-type systems against Table 1 of [Beverland et al., arXiv:2511.15177v1](https://arxiv.org/html/2511.15177v1). Both use `R = d` noisy cycles and one noise-free closing cycle.

| System | Noisy cycles | Detectors M | Logical observables K | Expanded fault weight N | Compressed columns N-tilde |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gross, BB(12)-circuit | 12 | 936 | 12 | 184,896 | 8,785 |
| Two-gross, BB(18)-circuit | 18 | 2,736 | 12 | 554,688 | 26,209 |

Here `N` is measured in units of `p/15`, not the number of gate locations. The compressed count includes one all-zero column. Matching these sizes is a necessary consistency check; it does not by itself prove that two fault models are identical.

## Detector-signature fingerprint

The expected two-gross Z-type fingerprint is:

```text
8002fb86c839b990f501ceb0237149c7a94ecc77b55b9f76da4cd052d8689769
```

The checksum covers 26,208 nonzero detector-signature entries and a zero-detector class of weight 20,736. It records detector responses and multiplicities; it does **not** include the logical-observable labels. Thus matching it is not a complete equality check of the detector-and-logical-action model.

To reproduce the fingerprint, follow `fingerprint()` in `verify_stim.py`:

1. Build the two-gross X-basis memory circuit with `R = 18` and `p = 10^-9`.
2. Obtain its undecomposed, flattened Stim detector error model.
3. For each error instruction with a nonempty detector response, sort its zero-based detector indices and format `d1,d2,...:u`, where `u` is the instruction's probability divided by `p/15`, rounded to the nearest integer. Logical targets are omitted. Detector indices are cycle-major, with X checks in upstream `decoder_setup.py` order.
4. Sort these strings lexicographically. Append `zero:u0`, where `u0 = 554688 - sum(u)`.
5. Join the lines with `\n`, with no final newline, encode as UTF-8, and take SHA-256.

For an independent comparison, derive the same detector responses and multiplicities from the upstream single-fault list with the same parameters and boundary. A full model comparison must also retain and compare logical actions.

## Interpretation

A silent, nontrivial 16-fault pattern can be partitioned into two eight-fault patterns. By linearity they have the same detector response and different logical actions. Any decoder using that response must fail on at least one of the two. This establishes a decoder-independent obstruction to correcting every eight-fault pattern. It does not establish that eight is the first failing weight: that additionally needs a lower-bound proof and an appropriate decoder.

The public witnesses tighten an upper bound. They do not contradict the statement `d_circ <= 18`, which remains true. Bravyi et al. explicitly noted that their bound might not be tight. Any application that assumes equality with 18 or a first failure at weight nine needs that assumption re-examined for the same circuit and noise convention. The change in practical failure probability cannot be inferred from distance alone.

The separately reported exhaustive search excluding all smaller logical errors is outside this release. No claim of reproducing that search, establishing novelty, or measuring performance at a given physical error rate is made by these scripts.

## References

- Bravyi et al., [*High-threshold and low-overhead fault-tolerant quantum memory*](https://arxiv.org/html/2308.07915v2), Nature 627, 778–782 (2024). Circuit construction and published upper bound, Table 1.
- Beverland et al., [*Fail fast: techniques to probe rare events in quantum error correction*](https://arxiv.org/html/2511.15177v1). System sizes in Table 1 and distance search in Table 2.
- Yoder et al., [*Tour de gross: A modular quantum computer based on bivariate bicycle codes*](https://arxiv.org/abs/2506.03094). Related architecture and schedule.
- Bravyi et al., [BivariateBicycleCodes source](https://github.com/sbravyi/BivariateBicycleCodes/tree/fa77e3333d3ec44c79d8f914dd24c040d1da471b).
