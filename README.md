# Secure State Estimation for Cyber-Physical Systems

This repository contains the MATLAB implementation for the secure state estimation of a dynamic Cyber-Physical System (CPS) subject to sparse, constant sensor attacks[cite: 2]. The project focuses on online state tracking and exact attack estimation using optimization-based sparse observers.

**System Model**
The unforced dynamics and output measurements of the CPS are modeled as follows:
* State equation: `x(k+1) = A * x(k)`
* Output equation: `y(k) = C * x(k) + a`

**Implemented Algorithms**
Because the system dynamics overlap with the constant attack, classical Luenberger observers fail[cite: 2]. To guarantee secure state estimation, the following sparse online optimization algorithms are implemented and compared:
* Online Proximal Gradient Descent (O-PGD)
* Online Proximal Alternating Minimization (O-PAM)
* Aggregated Online PGD (AO-PGD)
* Aggregated Online PAM (AO-PAM)

**Evaluated Metrics**
The performance and behavior of the algorithms are analyzed using four key metrics:
* **Support attack error:** Evaluates the discrepancy between the estimated and true attack support.
* **False positives:** Counts uncorrupted sensor channels mistakenly classified as attacked.
* **False negatives:** Counts active malicious injections missed by the estimator.
* **State estimation error:** Quantifies the normalized distance between the estimated state and the true state trajectory.

**Repository Contents**
* `CPS_part1.pdf`: The complete report detailing the theoretical analysis, dynamic scheduling, hyperparameter tuning, and graphical results.
* `main.m`: The primary MATLAB script that executes all four algorithms and generates the comparative plots.
* `dynamic_CPS_data.mat`: The dataset containing the system matrices (`A`, `C`), initial state (`x0`), and the true sparse attack vector (`a`) required to simulate the environment.

**How to Run**
1. Ensure `dynamic_CPS_data.mat` is located in the same directory as the scripts.
2. Run `main.m` in MATLAB.
3. The script will automatically compute the estimates and generate the performance graphs across all iterations.

**Authors**
* Buson Daniele
* Falco Abramo Calogero
* Franzoni Tommaso
* Martinelli Alex
Final Project for the Modeling and Control of Cyber-Physical Systems (CPSs) course at Politecnico di Torino.
