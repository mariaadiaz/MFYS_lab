Example data and core matlab functions to run the GUI used for the paper titled "A muscle synergy based method to improve robot-assisted movements"

# Muscle Synergy-Based Analysis GUI

## 1. **Overview of the Graphical User Interface (GUI)**

This tool allows users to compare muscle synergies to a set of reference muscle synergies, facilitating the analysis of muscle coordination patterns. It also supports optimization of one assistance parameter (e.g., force) to enhance human-robot interaction, enabling users to apply a gradient descent optimization strategy to determine the optimal assistance level for improved performance.

- **GUI Name**: Muscle Synergy-Based Analysis GUI
- **Version**: 1.0
- **Operating System**: Windows, macOS, or Linux
- **Developed in**: MATLAB (R2020b)

---

## 2. **How to Start Using the GUI**

Before using the GUI, ensure you have the reference muscle synergies prepared (denoted as `H_ref` and `W_ref`). This GUI processes datasets collected using the DELSYS system, with data stored in CSV format. Ensure your data file meets the following criteria:

- Row 4 contains the variable names.
- The file includes the sample frequency of the trial.

> _If your data structure differs, you may modify the code view in the .mlapp file accordingly._

### Example Data Files Provided:
- Reference muscle synergies in `.mat` format
- Sample data from a participant performing the task with three assistance levels (20, 100, and 200N)

**Input Requirements**:
- Two `.mat` files containing reference synergies (`W_ref` and `H_ref`)
- A CSV file with EMG data and sampling frequency

**Notes**:
- When **Both** is selected in outputs and **Calculate Muscle Synergies** is pressed, the GUI performs two calculations:
  - _Estimates_ `W_est` using `H_ref`
  - _Estimates_ `H_est` using `W_ref`
- For **Estimate H** or **Estimate W** individually, the GUI will output only for the selected option.

---

## 3. **User Guide (Step-by-Step)**

Follow these steps to load data, select outputs, and perform optimizations:

1. **Load the Reference Synergies**
   - Click **Load Reference Synergies**
   - Select both `.mat` files containing reference synergies simultaneously
   - A confirmation message will display: “Reference synergies were loaded”

2. **Select Desired Outputs**
   - Choose from **Estimate H**, **Estimate W**, or **Both**

3. **Use Optimization (Optional)**
   - Mark **Gradient Descent** checkbox
   - Input the force or parameter value for tuning; optimization will adjust this to improve performance

4. **Load a New Trial**
   - Click **Load New Trial**
   - A confirmation message “_filename_ was loaded” will appear if loaded successfully

5. **Calculate Muscle Synergies**
   - Press **Calculate Muscle Synergies** to start the estimation

**Additional Notes**:

- _For multiple trials_: 
  1. Load the new trial
  2. Select desired outputs
  3. Press **Calculate Muscle Synergies**

- _For gradient descent optimization_:
  1. Load a new trial
  2. Select desired outputs
  3. Adjust force (or parameter) as needed
  4. Press **Calculate Muscle Synergies**

After two iterations, a “_Next Value_” prompt will display a suggested parameter for further optimization to minimize the objective function.

---

## 4. **Features**

The GUI provides several metrics to evaluate the similarity between estimated and reference muscle synergies:

- **Total Abs Diff**:
  - Sum of absolute differences between `W_est` and `W_ref` for each synergy across all muscles
  - Individual muscle differences shown under **Diff (W est)** tab

- **Total ED (Euclidean Distance)**:
  - Point-by-point Euclidean distance between `H_est` and `H_ref`
  - Displayed under **ED (H est)** tab

- **Total SSV (Similarity of Synergy Vector)**:
  - Sum of highest similarity values (SSV) between `W_est` and `W_ref`, max value of 3
  - Displayed in color-coded format under **SSV (W est)**:
    - Dark red = high similarity (0.8–1)
    - Dark blue = low similarity (< 0.5)
