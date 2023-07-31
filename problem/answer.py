import sys
from typing import Any

sys.path.append("../")
from utils.challenge_2023 import ChallengeSampling

challenge_sampling = ChallengeSampling(noise=True)

"""
####################################
add codes here
####################################
"""

import time
import numpy as np
from openfermion.transforms import jordan_wigner
from openfermion.utils import load_operator
from quri_parts.algo.optimizer import Adam, OptimizerStatus
from quri_parts.core.estimator.gradient import parameter_shift_gradient_estimates
from quri_parts.core.measurement import bitwise_commuting_pauli_measurement
from quri_parts.core.sampling.shots_allocator import (
    create_equipartition_shots_allocator,
)
from quri_parts.core.state import ParametricCircuitQuantumState, ComputationalBasisState
from quri_parts.openfermion.operator import operator_from_openfermion_op
from utils.challenge_2023 import TimeExceededError
from quri_parts.circuit.circuit_parametric import UnboundParametricQuantumCircuit
from quri_parts.circuit import LinearMappedUnboundParametricQuantumCircuit
from quri_parts.circuit import QuantumCircuit
from qiskit import QuantumCircuit as QuantumCircuitQiskit
from qiskit.quantum_info import Operator
from quri_parts.qiskit.circuit import circuit_from_qiskit

def cost_fn(hamiltonian, parametric_state, param_values, estimator):
    estimate = estimator(hamiltonian, parametric_state, [param_values])
    return estimate[0].value.real

from qiskit import QuantumCircuit as QuantumCircuitQiskit
from qiskit.quantum_info import Operator
import numpy as np

def qiskit_gates(n_sites=4, trotter_steps=1, parameters=None, barrier=False, 
                  t=1.0, delta_t=1.0, U=3.0):      
    '''
    delta_t: time difference parameter.
    U: on-site Coulomb repulsion between pair of electrons occupying the same lattice site.
    t: tunneling amplitude parameter from the Hopping term.
    '''       
    n_qubits = 2 * n_sites
    qc = QuantumCircuitQiskit(n_qubits)  
    total_no_of_gates = int(n_sites + (n_qubits) // 2) #+ 4*(n_sites-1))
    num_parameters = trotter_steps * total_no_of_gates
    
    if parameters is None:
        parameters = np.random.rand(num_parameters)*2*np.pi*0.001
        
    theta = t * delta_t / 2
    phi = 1j * U * delta_t
    
    def XX_gate(parameter):
        return Operator([[np.cos(parameter*theta), 0, 0, 1j * np.sin(parameter*theta)],
                         [0, np.cos(parameter*theta), 1j * np.sin(parameter*theta), 0],
                         [0, 1j * np.sin(parameter*theta), np.cos(parameter*theta), 0],
                         [1j * np.sin(parameter*theta), 0, 0, np.cos(parameter*theta)]])

    def YY_gate(parameter):
        return Operator([[np.cos(parameter*theta), 0, 0, -1j * np.sin(parameter*theta)],
                         [0, np.cos(parameter*theta), 1j * np.sin(parameter*theta), 0],
                         [0, 1j * np.sin(parameter*theta), np.cos(parameter*theta), 0],
                         [-1j * np.sin(parameter*theta), 0, 0, np.cos(parameter*theta)]])

    def Givens_Rotations(parameter):
        return Operator([[1, 0, 0, 0],
                         [0, np.cos(parameter/2), -np.sin(parameter/2), 0],
                         [0, np.sin(parameter/2), np.cos(parameter/2), 0],
                         [0, 0, 0, 1]]) 
    
    def on_site_gate(parameter, i, j):
        unitary_matrix = Operator([[1, 0, 0, 0],
                                   [0, 1, 0, 0],
                                   [0, 0, 1, 0],
                                   [0, 0, 0, np.exp(parameter*phi)]])
        qc.unitary(unitary_matrix, [i, j], label=f"Onsite θ[{parameter_index}]")
    
    def hopping_gate(parameter, i, j):
        qc.unitary(YY_gate(parameter), [i, j], label=f"YY Hop. θ[{parameter_index}]")
        qc.unitary(XX_gate(parameter), [i, j], label=f"XX Hop θ[{parameter_index}]")
        
    def hopping_gate2(parameter, i, j):
        unitary_matrix = np.array([[1, 0, 0, 0],
                                   [0, np.cos(parameter/2), -1j*np.sin(parameter/2), 0],
                                   [0, -1j*np.sin(parameter/2), np.cos(parameter/2), 0],
                                   [0, 0, 0, 1]])
        qc.unitary(unitary_matrix, [i, j], label=f"Hopping θ[{parameter_index}]")

    def get_angle():
        nonlocal parameter_index
        angle = parameters[parameter_index]
        parameter_index += 1
        return angle

    parameter_index = 0
    for _ in range(trotter_steps):
        '''
        for i in range(1, n_qubits-1, n_sites):
            qc.unitary(Givens_Rotations(get_angle()), [i, i+1], label=f"Givens-Rotations θ[{parameter_index}]")

        for i in range(0, n_qubits-1, 2):
            qc.unitary(Givens_Rotations(get_angle()), [i, i+1], label=f"Givens-Rotations θ[{parameter_index}]")

        for i in range(1, n_qubits-1, n_sites):
            qc.unitary(Givens_Rotations(get_angle()), [i, i+1], label=f"Givens-Rotations θ[{parameter_index}]")
        '''    
        if barrier:
            qc.barrier()
        
        for i in range(n_sites):
            on_site_gate(get_angle(), i, i + n_sites)
            
        if barrier:
            qc.barrier()
            
        #for i in range(0, n_sites-1):
            #hopping_gate(get_angle(), i, i+1)
            #hopping_gate(get_angle(), i+n_sites, i+1+n_sites)
            
        for i in range(0, n_qubits, 2):
            hopping_gate(get_angle(), i, i+1)
            #hopping_gate2(get_angle(), i, i + 1)
            
        if barrier:
            qc.barrier()
    
    return qc    

def quri_gates(n_sites, trotter_steps, params, option, t=1.0, delta_t=1.0, U=3.0):
    n_qubits = 2*n_sites
    if option == 1:
        circuit = UnboundParametricQuantumCircuit(n_qubits)
    elif option == 2:
        circuit = QuantumCircuit(n_qubits)
    theta = t * delta_t / 2
    phi = 1j * U * delta_t
    index = 0
    
    sqrt_iS = np.array([[1, 0, 0, 0],
                        [0, 1/np.sqrt(2), 1j/np.sqrt(2), 0],
                        [0, 1j/np.sqrt(2), 1/np.sqrt(2), 0],
                        [0, 0, 0, 1]])
        
    def XX_gate(parameter, i, j):
        circuit.add_H_gate(i)
        circuit.add_H_gate(j)
        circuit.add_CNOT_gate(i, j)
        if option == 1:
            circuit.add_ParametricRZ_gate(j)
        elif option == 2:
            circuit.add_RZ_gate(index=j, angle=theta*parameter) 
        circuit.add_CNOT_gate(i, j)
        circuit.add_H_gate(i)
        circuit.add_H_gate(j)

    def YY_gate(parameter, i, j):
        circuit.add_Sdag_gate(i)
        circuit.add_Sdag_gate(j)
        circuit.add_H_gate(i)
        circuit.add_H_gate(j)
        circuit.add_CNOT_gate(i, j)
        if option == 1:
            circuit.add_ParametricRZ_gate(j)
        elif option == 2:
            circuit.add_RZ_gate(index=j, angle=theta*parameter) 
        circuit.add_CNOT_gate(i, j)
        circuit.add_H_gate(i)
        circuit.add_H_gate(j)
        circuit.add_S_gate(i)
        circuit.add_S_gate(j)

    '''
    def on_site_gate(i, j):
        circuit.add_ParametricRZ_gate(i)
        circuit.add_ParametricRZ_gate(j)
        circuit.add_ParametricRX_gate(i) 
        circuit.add_RX_gate(index=j, angle=-np.pi/2) 
        circuit.add_Z_gate(i)
        circuit.add_UnitaryMatrix_gate([i,j], sqrt_iS)
        circuit.add_Z_gate(i)
        circuit.add_ParametricRX_gate(j)
        circuit.add_UnitaryMatrix_gate([i,j], sqrt_iS)
        circuit.add_ParametricRX_gate(i)
        circuit.add_RX_gate(index=j, angle=np.pi/2)
    '''
    
    def onsite_gate(parameter, i, j):
        unitary_matrix = np.array([[1, 0, 0, 0],
                                   [0, 1, 0, 0],
                                   [0, 0, 1, 0],
                                   [0, 0, 0, np.exp(phi*parameter)]])
        circuit.add_UnitaryMatrix_gate([i, j], unitary_matrix)
    
    def hopping_gate(parameter, i, j):
        if option == 1:
            YY_gate(parameter, i, j)
            XX_gate(parameter, i, j)
        elif option == 2:
            unitary_matrix = np.array([[1, 0, 0, 0],
                                       [0, np.cos(parameter/2), -1j*np.sin(parameter/2), 0],
                                       [0, -1j*np.sin(parameter/2), np.cos(parameter/2), 0],
                                       [0, 0, 0, 1]])
            circuit.add_UnitaryMatrix_gate([i, j], unitary_matrix)
               
    for _ in range(trotter_steps):
        
        for i in range(n_sites):
            #on_site_gate(i, i + n_sites)
            onsite_gate(params[index], i, i + n_sites)
            index+=1
            
        for i in range(0, n_qubits, 2):
            hopping_gate(params[index], i, i + 1)
            index+=1
    
    return circuit

def rotosolve(cost, params):
    '''This function defines the rotosolve optimizer.'''
    #phi: float = float(np.random.uniform(0, np.pi))
    phi, k = 0, 0
    for d in range(len(params)):
        params[d] = phi
        M_0 = cost(params)
        params[d] = phi + np.pi / 2.0
        M_0_plus = cost(params)
        params[d] = phi - np.pi / 2.0
        M_0_minus = cost(params)
        a = np.arctan2(2.0 * M_0 - M_0_plus - M_0_minus, M_0_plus - M_0_minus)  
        params[d] = phi - np.pi / 2.0 - a + 2*np.pi*k
        if params[d] <= -np.pi:
            params[d] += 2 * np.pi
    return params  

def vqe(hamiltonian, parametric_state, estimator, params, optimizer):
    energies = []
    print('Computing VQE...\n')
    
    def c_fn(params):
        return cost_fn(hamiltonian, parametric_state, params, estimator)
    
    #for iteration in range(1, 50+1, 1):
    iteration = 0
    while True:
        try:
            start_time = time.time()
            params = optimizer(c_fn, params)
            _energy = c_fn(params)
            print(f"iteration {iteration+1}")
            print(_energy)
            end_time = time.time()
            execution_time = end_time - start_time
            print("Execution Time:", execution_time, "seconds")
            energies.append(_energy)
            iteration+=1
        except TimeExceededError as e:
            print(str(e))
            return _energy
    
    return _energy, iteration, energies
    
class RunAlgorithm:
    def __init__(self) -> None:
        challenge_sampling.reset()

    def result_for_evaluation(self) -> tuple[Any, float]:
        energy_final = self.get_result()
        qc_time_final = challenge_sampling.total_quantum_circuit_time

        return energy_final

    def get_result(self) -> float:
        """
        ####################################
        add codes here
        ####################################
        """
        n_sites, trotter_steps = 4, 1
        n_qubits = 2*n_sites
        ham = load_operator(
            file_name=f"{n_qubits}_qubits_H",
            data_directory="../hamiltonian",
            plain_text=False,
        )
        # Jordan-Wigner Transformation:
        jw_hamiltonian = jordan_wigner(ham)
        # Spin qubit-equivalent Hamiltonian:
        hamiltonian = operator_from_openfermion_op(jw_hamiltonian)
        
        # Parameter Initialization:
        total_no_of_gates = int(n_sites + (n_qubits) // 2) #+ 4*(n_sites-1))
        init_params = np.random.rand(trotter_steps * total_no_of_gates)*2*np.pi*0.001
        
        ##### Circuit Ansatz: #####
        
        hf_gates = ComputationalBasisState(n_qubits, bits=0b00001111).circuit.gates
        # Option1 - Quri & UnboundParametricQuantumCircuit
        #hf_circuit = quri_gates(n_sites, trotter_steps, init_params, option = 1)
        #hf_circuit.combine(hf_gates)
        
        # Option2 - Quri & LinearMappedUnbound...
        #hf_circuit = LinearMappedUnboundParametricQuantumCircuit(n_qubits).combine(hf_gates)
        #hw_ansatz = quri_gates(n_sites, trotter_steps, init_params, option=2)
        #hf_circuit.combine(hw_ansatz)
        
        # Option3 - Qiskit & LinearMappedUnbound...
        hf_circuit = LinearMappedUnboundParametricQuantumCircuit(n_qubits).combine(hf_gates)
        hw_ansatz = qiskit_gates(n_sites, trotter_steps, init_params)
        hf_circuit.combine(circuit_from_qiskit(qiskit_circuit = hw_ansatz).gates)
        
        ##### Circuit Ansatz: #####
        
        # Parametric state:
        parametric_state = ParametricCircuitQuantumState(n_qubits, hf_circuit)

        hardware_type = "sc"
        shots_allocator = create_equipartition_shots_allocator()
        measurement_factory = bitwise_commuting_pauli_measurement
        n_shots = 10**3
        
        # Sampling estimator to estimate the expectation value of the Hamiltonian operator:
        sampling_estimator = (
            challenge_sampling.create_concurrent_parametric_sampling_estimator(
                n_shots, measurement_factory, shots_allocator, hardware_type
            )
        ) 
        # VQE Simulation:
        energy_final = vqe(
            hamiltonian,
            parametric_state,
            sampling_estimator,
            init_params,
            rotosolve,
        )
        print(f"iteration used: {iteration}")
        return energy_final

if __name__ == "__main__":
    run_algorithm = RunAlgorithm()
    print(run_algorithm.get_result())
