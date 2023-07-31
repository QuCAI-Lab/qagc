<div align="center">
  <a href="https://qucai-lab.github.io/"><img src="https://github.com/QuCAI-Lab/qucai-lab.github.io/blob/main/assets/QuCAI-Lab.png" height="500" width="500" /></a>
</div>

<div align="center">
  <h1> <a href="https://github.com/QunaSys/quantum-algorithm-grand-challenge/tree/main"> Quantum Algorithm Grand Challenge</a></h1>
  <h2> Applying the HVA method to find the Fermi-Hubbard Ground State </h2>
</div>
<br>

<div align="center">
  <b>Author: <a target="_blank" href="https://github.com/camponogaraviera">Lucas Camponogara Viera</a></b>
</div>

# Presentation

Refer to the [Presentation.ipynb](Presentation.ipynb) notebook for a theoretical background and derivation of the Hopping gate and Onsite gate using the Jordan–Wigner mapping and the Suzuki-Trotter formula. Gradient-free and gradient-based optimizations were compared and the circuit ansatz was implemented using both [Qiskit](https://qiskit.org/) and [QURI Parts](https://quri-parts.qunasys.com/).

# Project Description

Various approaches have been proposed to create effective hardware-efficient and chemistry-inspired ansatz. One such ansatz is the trotterized UCCSD ansatz which, although based on the accurate Coupled-cluster theory, suffers from inconsistency under low-order trotterization steps [1] when shallow circuits friendly to NISQ devices provide a poor approximation to the UCCSD wavefunction. 

On the optimization side, gradient-based algorithms for the optimization of variational quantum circuits have also been explored, such as the Bayesian model gradient descent (BayesMGD) [2], the model gradient descent (MGD) [3], the [SPSA](https://quri-parts.qunasys.com/quri_parts/algo/quri_parts.algo.optimizer.html#quri_parts.algo.optimizer.SPSA), [Adam](https://quri-parts.qunasys.com/quri_parts/algo/quri_parts.algo.optimizer.html#quri_parts.algo.optimizer.Adam), [L-BFGS](https://quri-parts.qunasys.com/quri_parts/algo/quri_parts.algo.optimizer.html#quri_parts.algo.optimizer.LBFGS), [NFT](https://quri-parts.qunasys.com/quri_parts/algo/quri_parts.algo.optimizer.html#quri_parts.algo.optimizer.NFT), and the ADAPT-VQE [4] which strongly relies on the parameter-shift rule to compute the gradient of the gates in order to optimize the Euler angles and change the circuit structure on-the-fly [5]. Gradient-free optimization algorithms include the Rotosolve, Rotoselect [6], and the NelderMead method, to name a few. The mainstream literature states that gradient-based methods are more efficient with faster convergence if gradients can be computed directly [7] [8]. 

The characteristic of the Hubbard model suggests the application of the Hamiltonian variational method [9] which uses terms of the Hamiltonian to propose the circuit ansatz with a smaller circuit depth than the unitary coupled cluster method (UCC). 

In this solution, gradient-free and gradient-based optimizations were compared, and the final algorithm was chosen based on its computational time (speed) and convergence.

\[1] Grimsley, H. R.; Claudino, D.; Economou, S. E.; Barnes, E.; Mayhall, N. J. Is the trotterized uccsd ansatz chemically well-defined? [J. Chem. Theory Comput. 2020, 16, 1](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b01083).

\[2] Stanisic, S., Bosse, J.L., Gambetta, F.M. et al. Observing ground-state properties of the Fermi-Hubbard model using a scalable algorithm on a quantum computer. [Nat Commun 13, 5743 (2022).](https://www.nature.com/articles/s41467-022-33335-4)

\[3] Sung, K. J. et al. Using models to improve optimizers for variational quantum algorithms. Quantum Sci. Technol. 5, 044008 (2020).

\[4] Harper R. Grimsley, Sophia E. Economou, Edwin Barnes, Nicholas J. Mayhall, "An adaptive variational algorithm for exact molecular simulations on a quantum computer". [Nat. Commun. 2019, 10, 3007](https://www.nature.com/articles/s41467-019-10988-2).

\[5] Lucas Camponogara Viera and José Paulo Marchezi. "Extending Adaptive Methods for Finding an Optimal Circuit Ansatze in VQE Optimization". [QHack2022 Open Hackathon Project - Quantum Chemistry Challenge.](https://github.com/QuCAI-Lab/qhack2022-hackeinberg-project)

\[6] Mateusz Ostaszewski, Edward Grant, Marcello Bendetti. "Structure optimization for parameterized quantum circuits." [Quantum 5, 391 (2021).](https://quantum-journal.org/papers/q-2021-01-28-391/)

\[7] J. Li, X. Yang, X. Peng, and C.-P. Sun, [Phys. Rev. Lett. 118, 150503 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.150503).

\[8] K. Mitarai, M. Negoro, M. Kitagawa, and K. Fujii, [Phys. Rev. A 98, 032309 (2018)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.032309).

\[9] Wecker, D., Hastings, M. B. & Troyer, M. Progress towards practical quantum variational algorithms. [Phys. Rev. A 92](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.92.042303), 042303 (2015).

\[10] Chris Cade, Lana Mineh, Ashley Montanaro, and Stasja Stanisic. Strategies for solving the Fermi-Hubbard model on near-term quantum computers. [Phys. Rev. B 102, 235122](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.235122).

\[11] Jiang, Z., Sung, K. J., Kechedzhi, K., Smelyanskiy, V. N. & Boixo, S. Quantum algorithms to simulate many-body physics of correlated fermions. Phys. Rev. Appl. 9, 044036 (2018).
 
# Environment

```bash
conda create -yn qagc-cpu python==3.10.9 && conda activate qagc-cpu
```

# Dependencies

```bash
python -m pip install -r requirements.txt
```

# Quick run

```bash
cd problem && python answer.py
```
