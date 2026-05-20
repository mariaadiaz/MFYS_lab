# HILO / BILO - Real-time EMG Synergy and Bayesian Optimization

This project is built on top of the **Delsys Trigno API** (Python SDK provided by Delsys, 2023). The original vendor sample has been extended with a processing and optimization layer that runs in real time during a human-in-the-loop (HILO / BILO) gait study.

The Delsys API handles sensor pairing, streaming, and the .NET bridge (see [AeroPy/](AeroPy/) and [resources/](resources/)). Everything in [Functions/](Functions/), [DataCollector/](DataCollector/), and [Plotter/](Plotter/) was added on top.

> More information about the Delsys API itself can be found in [README delsys.md](README%20delsys.md).

---

## `Functions/synergies.py` — brief description

The file contains a single class, `Processingcollector`, with the following methods:

- **`process_emg(df, fs, max_*)`** — Builds the linear envelope of each muscle EMG signal, normalizes by the running per-muscle maximum, slices the last 20 s of data into gait cycles using gastrocnemius peaks, averages cycles, and runs NMF (4 synergies) to obtain `W` (synergy weights) and `H` (activation patterns). Also returns muscle power (Welch PSD) and VAF.

- **`organize_synergies(W_norm, W2_norm, H2, n)`** — Reorders the synergies of a second window so that each column best matches the reference window, using cosine similarity. Required before comparing the two windows.

- **`compare_muscle_outputs(power1, power2, W1, W2, H1, H2, n)`** — Compares synergies and activations between the min-2 and min-3 windows of the same stride-frequency iteration. Outputs mean SSV (angle between matched synergies), mean Euclidean distance between activations, mean power, and the predicted energy cost (EC).

- **`estimation(mean_ssv, mean_dist, mean_power)`** — Linear regression that maps the three motor-variability features into a predicted EC value. Coefficients are taken from Díaz et al., 2024.

- **`bayesian_optimization(sf_vals, y_vals)`** — Gaussian-process regression (RBF kernel) over previously tested stride frequencies, combined with an expected-improvement acquisition function, to suggest the next SF to test (search space 75–126 % of preferred SF).

- **`expected_improvement(x, gp_model, best_y)`** — Standard EI formula used by the Bayesian optimizer.

- **`Partial_results_callback` / `results_callback`** — Helper callbacks that expose the current results to the UI controller.

---

## Important notes before running

### 1. Update the output file names
At the end of the session the controller saves the data to hardcoded files (currently `Testing.pkl` and `Testing.txt` at the project root, see [DataCollector/CollectDataController.py](DataCollector/CollectDataController.py)). These names **must be changed for every new participant / session**, otherwise previous results will be overwritten.

### 2. The code expects exactly 8 Delsys sensors
Sensors must be connected and placed on the following muscles (the pair number selects the muscle label):

| Pair number | Muscle |
|---|---|
| 1 / 9 | Right Rectus Femoris |
| 2 / 10 | Right Gastrocnemius Lateralis |
| 3 / 11 | Right Gastrocnemius Medialis |
| 4 / 12 | Right Tibialis Anterior |
| 5 / 13 | Left Rectus Femoris |
| 6 / 14 | Left Gastrocnemius Lateralis |
| 7 / 15 | Left Gastrocnemius Medialis |
| 8 / 16 | Left Tibialis Anterior |

Removing, adding, or reassigning a sensor will break the processing pipeline, because muscle columns are selected by name and the synergy analysis assumes all eight channels.

### 3. Offline mode still requires the Delsys system
It is possible to re-run the pipeline on previously collected data (`.pkl` files), but the Delsys base station **must still be connected via USB** and the API must be able to load successfully. This is because the application initializes the Delsys .NET bridge at startup, regardless of whether the data source is live or recorded.
