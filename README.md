# Secure State Estimation for Cyber-Physical Systems

This repository contains the MATLAB implementation for the secure state estimation of a dynamic Cyber-Physical System (CPS) subject to sparse, constant sensor attacks[cite: 2]. The project focuses on online state tracking and exact attack estimation using optimization-based sparse observers[cite: 2].

**System Model**
The unforced dynamics and output measurements of the CPS are modeled as follows[cite: 2]:
* State equation: `x(k+1) = A * x(k)`[cite: 2]
* Output equation: `y(k) = C * x(k) + a`[cite: 2]

**Implemented Algorithms**
Because the system dynamics overlap with the constant attack, classical Luenberger observers fail[cite: 2]. To guarantee secure state estimation, the following sparse online optimization algorithms are implemented and compared[cite: 2]:
* Online Proximal Gradient Descent (O-PGD)[cite: 2]
* Online Proximal Alternating Minimization (O-PAM)[cite: 2]
* Aggregated Online PGD (AO-PGD)[cite: 2]
* Aggregated Online PAM (AO-PAM)[cite: 2]

**Evaluated Metrics**
The performance and behavior of the algorithms are analyzed using four key metrics[cite: 2]:
* **Support attack error:** Evaluates the discrepancy between the estimated and true attack support[cite: 2].
* **False positives:** Counts uncorrupted sensor channels mistakenly classified as attacked[cite: 2].
* **False negatives:** Counts active malicious injections missed by the estimator[cite: 2].
* **State estimation error:** Quantifies the normalized distance between the estimated state and the true state trajectory[cite: 2].

**Repository Contents**
* `CPS_part1.pdf`: The complete report detailing the theoretical analysis, dynamic scheduling, hyperparameter tuning, and graphical results[cite: 2].
* `main.m`: The primary MATLAB script that executes all four algorithms and generates the comparative plots.
* `dynamic_CPS_data.mat`: The dataset containing the system matrices (`A`, `C`), initial state (`x0`), and the true sparse attack vector (`a`) required to simulate the environment[cite: 2].

**How to Run**
1. Ensure `dynamic_CPS_data.mat` is located in the same directory as the scripts.
2. Run `main.m` in MATLAB.
3. The script will automatically compute the estimates and generate the performance graphs across all iterations.
