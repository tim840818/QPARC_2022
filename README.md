# Second-Quantization Comparison of $H_2$ Molecule
 
In this work, we prepare two kinds of Second-Quantized Hamiltonian from `PySCF` and Qulacs for comparison. The prepared Hamiltonians are then transfered to JW-Hamiltonians from `openfermion`, and the ground-state energy are calculated by implementing the Hamiltonians into quantum circuits constructed by `Qulacs`. The calculations are done by VQE algorithm, and the quantum circuits are constructed with hardware efficient ansatz shown in [Qulacs](http://docs.qulacs.org/en/latest/apply/6.2_vqe.html). The results are compared with SCF results from `PySCF` and FCI results from `Qulacs`.

For the results, we get the following ground state energies:\\
 Calculation method | Time  
:-------------------|:-----:
PySCF Hamiltonian   | -1.6820 
Qulacs Hamiltonian  | -1.1178 
Classical SCF       | -1.0985 
Classical FCI       | -1.1254 

The ground state energy from `PySCF` Hamiltonian is about 0.6 Hartree lower than other results. There might be some parameters missing that leads to this difference.


## Packages

- `PySCF`: creating 2nd-quantized Hamitonians and performing classical scf calculations.
- `openfermion`: transforming 2nd-quantized Hamiltonians into JW-Hamiltonians.
- `Qulacs`: creating 2nd-quantized Hamitonians, performing fci calculations, and performing quantum simulations with quantum circuits.

## Presentation

The sub-folder `presentation` contains:\\
- `Second Quantization Study of H2.pptx` as our powerpoint presentation.
- `H2_qulacspyscf.ipynb` as our jupyter-notebook presentation.

## Code

The sub-folder `code` contains `H2_qulacspyscf.py` created from `H2_qulacspyscf.ipynb`.
