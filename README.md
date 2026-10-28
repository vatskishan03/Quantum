# Two-gross circuit-distance witnesses

**Sixteen gate faults can change stored quantum information without triggering any detector in the two-gross memory circuit studied here.** This repository provides the exact faults and runnable checks for both Z-type and X-type errors.

The result is **circuit distance at most 16**, tightening the published upper bound of 18 for the specified schedule. This concerns the measurement circuit; it does not change the code's `[[288,12,18]]` parameters. See [Bravyi et al., Table 1](https://arxiv.org/html/2308.07915v2) and the [verification details](VERIFICATION.md).

## Check it yourself

You need Git and Python 3.12 or newer. No quantum computer, GPU, account, or search software is needed.

On macOS or Linux, run:

```sh
git clone https://github.com/vatskishan03/Quantum.git
cd Quantum
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py
```

On Windows, create the environment with `py -m venv .venv` and use `.venv\Scripts\python.exe` in place of `python` for the last two commands. No activation is needed.

After setup, the entire check is one command:

```sh
python verify.py
```

A successful run ends with:

```text
VERIFIED: all four 16-fault witnesses pass (circuit distance <= 16).
```

The command checks all seven file hashes, runs both verification methods, and checks the expected detector-signature fingerprint. It exits with a nonzero status if any check fails. Once the dependencies are installed, verification runs offline.

Tested with Python 3.12.14 and 3.15.0, NumPy 2.5.4, and Stim 1.16.0. The same command runs automatically on GitHub for pushes and pull requests.

## What is checked

| Check | What it establishes |
| --- | --- |
| Stim simulation | Each of four witnesses uses 16 allowed CNOT faults, fires zero detectors, and changes a logical observable. |
| Separate NumPy simulation | The two original-schedule witnesses remain undetected with either one or two noise-free closing cycles. Their final data errors are nontrivial logical operators. |
| Subset enumeration | No nonempty proper subset of any supplied witness is undetected. |
| Model consistency | The gross and two-gross Z-type systems reproduce the sizes in [Beverland et al., Table 1](https://arxiv.org/html/2511.15177v1). |

The four witnesses cover both error types in the original schedule and in a variant that exchanges the roles of `A2` and `A3`.

**Scope:** These checks prove an upper bound of 16. Checking subsets of one witness does not rule out other, smaller witnesses. The separately reported exhaustive search is not included, so this repository does not reproduce the lower bound needed to establish distance exactly 16. It also does not measure logical failure rates at a particular physical error rate.

## Files

| File | Purpose |
| --- | --- |
| [verify_stim.py](verify_stim.py) | Stim checks for all four witnesses and model sizes. |
| [gen_bb_memory.py](gen_bb_memory.py) | Builds the memory circuits. |
| [verify_transcription.py](verify_transcription.py) | Separate NumPy check of the original schedule. |
| [witness_Ztype_16faults.json](witness_Ztype_16faults.json), [witness_Xtype_16faults.json](witness_Xtype_16faults.json) | Original-schedule faults. |
| [witness_A23swap_Ztype_16faults.json](witness_A23swap_Ztype_16faults.json), [witness_A23swap_Xtype_16faults.json](witness_A23swap_Xtype_16faults.json) | Variant-schedule faults. |

[VERIFICATION.md](VERIFICATION.md) documents the circuit, noise convention, fingerprint, expected results, and interpretation. [SHA256SUMS](SHA256SUMS) records the seven source-file hashes. [verify.py](verify.py) runs the full check.

## Attribution and license

The circuit construction and transcription follow [Bravyi et al.'s BivariateBicycleCodes implementation](https://github.com/sbravyi/BivariateBicycleCodes/tree/fa77e3333d3ec44c79d8f914dd24c040d1da471b). See [NOTICE](NOTICE) for the source and adaptations. Released under the [Apache License 2.0](LICENSE).
