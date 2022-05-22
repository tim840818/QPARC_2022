import numpy as np
import matplotlib.pyplot as plt

from openfermion.transforms import get_fermion_operator, jordan_wigner
from openfermion.linalg import get_sparse_operator
from openfermion.chem import MolecularData
from openfermionpyscf import run_pyscf
from scipy.optimize import minimize
from pyscf import fci, lo, gto, dft, scf, cc, df, ao2mo
from openfermion.ops import InteractionOperator

import qulacs
from qulacs import Observable
from qulacs.observable import create_observable_from_openfermion_text

from qulacs import QuantumState, QuantumCircuit
from qulacs.gate import CZ, RY, RZ, merge


def round2(A, d=1):
    return np.round(A, d)

def he_ansatz_circuit(n_qubit, depth, theta_list):
    """he_ansatz_circuit
    Returns hardware efficient ansatz circuit.

    Args:
        n_qubit (:class:`int`):
            the number of qubit used (equivalent to the number of fermionic modes)
        depth (:class:`int`):
            depth of the circuit.
        theta_list (:class:`numpy.ndarray`):
            rotation angles.
    Returns:
        :class:`qulacs.QuantumCircuit`
    """
    circuit = QuantumCircuit(n_qubit)
    for d in range(depth):
        for i in range(n_qubit):
            circuit.add_gate(merge(RY(i, theta_list[2*i+2*n_qubit*d]), RZ(i, theta_list[2*i+1+2*n_qubit*d])))
        for i in range(n_qubit//2):
            circuit.add_gate(CZ(2*i, 2*i+1))
        for i in range(n_qubit//2-1):
            circuit.add_gate(CZ(2*i+1, 2*i+2))
    for i in range(n_qubit):
        circuit.add_gate(merge(RY(i, theta_list[2*i+2*n_qubit*depth]), RZ(i, theta_list[2*i+1+2*n_qubit*depth])))

    return circuit

def cost(theta_list):
    state = QuantumState(n_qubit) #Prepare |00000>
    circuit = he_ansatz_circuit(n_qubit, depth, theta_list) #Construct quantum circuit
    circuit.update_quantum_state(state) #Operate quantum circuit on state
    return qulacs_hamiltonian.get_expectation_value(state) #Calculate expectation value of Hamiltonian

if __name__ == "__main__":
## construct 2nd-quantized Hamiltonian from PySCF
    mole_H2 = gto.M(
        atom = '''H  0 0 0; H 0 0 0.87''',
        basis = 'sto-3g', # '6-31g'
        # symmetry = True,
        # charge = 0,
    )
    H2_ao = mole_H2.ao_labels(); print(H2_ao)
    n_qubit = len(H2_ao) #; print(n_qubit)
## calculate 1-body/2-body ao/mo integrals 
    mocoeff = mole_H2.RHF().run().mo_coeff
    kin = mole_H2.intor('int1e_kin')
    nuc = mole_H2.intor('int1e_nuc')
    # ovlp = mole_H2.intor('int1e_ovlp')
    twobody_ao = mole_H2.intor('int2e')
    core_ao = kin + nuc
    core_mo = np.einsum('pi,pq,qj->ij', mocoeff, core_ao, mocoeff)
    twobody_mo = ao2mo.incore.full(twobody_ao, mocoeff)
## transform to JW
    mole_H2_Hamiltonian = InteractionOperator(0, core_mo, twobody_mo)
    fermionic_hamiltonian = get_fermion_operator(mole_H2_Hamiltonian)
    jw_hamiltonian = jordan_wigner(fermionic_hamiltonian)
## calculate groundstate energy (PySCF) and plot
    depth = n_qubit

    qulacs_hamiltonian = create_observable_from_openfermion_text(str(jw_hamiltonian))
    cost_history_pyscf = []
    init_theta_list = np.random.random(2*n_qubit*(depth+1))*1e-1
    cost_history_pyscf.append(cost(init_theta_list))
    method = "BFGS"
    options = {"disp": True, "maxiter": 40, "gtol": 1e-6}
    opt = minimize(cost, init_theta_list,
                method=method,
                callback=lambda x: cost_history_pyscf.append(cost(x)))
    plt.rcParams["font.size"] = 18
    plt.plot(cost_history_pyscf, color="red", label="PySCF-VQE")
    plt.xlabel("Iteration")
    plt.ylabel("Energy expectation value")
    plt.legend()
    plt.show()
## ------------------------------ ##
## construct 2nd-quantized Hamiltonian from Qulacs
    basis = "sto-3g"
    multiplicity = 1
    charge = 0
    distance  = 0.87
    geometry = [["H", [0,0,0]],["H", [0,0,distance]]]
    description  = "tmp"
    molecule = MolecularData(geometry, basis, multiplicity, charge, description)
    molecule = run_pyscf(molecule,run_scf=1,run_fci=1)
    n_qubit = molecule.n_qubits
    n_electron = molecule.n_electrons

    mole_H = molecule.get_molecular_hamiltonian()
## transform to JW
    fermionic_hamiltonian = get_fermion_operator(mole_H)
    jw_hamiltonian = jordan_wigner(fermionic_hamiltonian)
## calculate groundstate energy (Qulacs) and plot
    depth = n_qubit

    qulacs_hamiltonian = create_observable_from_openfermion_text(str(jw_hamiltonian))
    cost_history_qulac = []
    init_theta_list = np.random.random(2*n_qubit*(depth+1))*1e-1
    cost_history_qulac.append(cost(init_theta_list))
    method = "BFGS"
    options = {"disp": True, "maxiter": 40, "gtol": 1e-6}
    opt = minimize(cost, init_theta_list,
                method=method,
                callback=lambda x: cost_history_qulac.append(cost(x)))
    plt.rcParams["font.size"] = 18
    plt.plot(cost_history_qulac, color="blue", label="Qulacs-VQE")
    plt.xlabel("Iteration")
    plt.ylabel("Energy expectation value")
    plt.legend()
    plt.show()
## ------------------------------ ##
## plot results
    E_SCF = -1.0985 

    iters = max(len(cost_history_qulac), len(cost_history_pyscf))

    plt.plot(cost_history_pyscf, color="red", label="PySCF-VQE")
    plt.plot(cost_history_qulac, color="blue", label="Qulacs-VQE")
    plt.plot(range(iters), [E_SCF]*iters, linestyle="dashed", color="green", label="SCF")
    plt.plot(range(iters), [molecule.fci_energy]*iters, linestyle="dashed", color="black", label="FCI")

    plt.title("Ground state energy of H2 molecule", fontsize=16)
    plt.xlabel("Iteration", fontsize=12)
    plt.ylabel("Energy expectation value (Hartree)", fontsize=12)
    plt.legend(fontsize=13)
    plt.show()
