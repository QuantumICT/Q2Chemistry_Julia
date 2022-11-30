module QuantumCircuits

using SparseArrays, LinearAlgebra, TensorOperations

# gates
export nqubits, positions, mat, ordered_positions, ordered_mat, change_positions, shift, differentiate
export parameters, nparameters, active_parameters, activate_parameter!, activate_parameters!, deactivate_parameter!, deactivate_parameters!, reset_parameters!
export Gate, ParametricGate, QuantumGate, AdjointQuantumGate, gate, XGate, YGate, ZGate, SGate, HGate, TGate, sqrtXGate, sqrtYGate
export SWAPGate, iSWAPGate, CZGate, CNOTGate, CONTROLGate
export TOFFOLIGate, FREDKINGate, CONTROLCONTROLGate

# parametric gates
export RxGate, RyGate, RzGate, PHASEGate
export CRxGate, CRyGate, CRzGate, CPHASEGate, FSIMGate
export CCPHASEGate


# quantum channel
export AbstractQuantumMap, QuantumMap, ordered_supermat, kraus_matrices, is_tp
export AmplitudeDamping, PhaseDamping, Depolarizing


# circuit 
export QMeasure, QSelect, QCircuit

# hamiltonian
export QubitsTerm, oplist, coeff, QubitsOperator, matrix, simplify

abstract type QuantumOperation end

# auxiliary
include("auxiliary/tensorops.jl")


# elementary gate matrices
include("elemops.jl")

using QuantumCircuits.Gates


# gate operations
include("gates/gates.jl")
include("gates/generic.jl")
include("gates/parametric_gates.jl")
include("gates/onebody.jl")
include("gates/twobody.jl")
include("gates/threebody.jl")
include("gates/gatediff.jl")

# quantum channel
include("channels/channels.jl")
include("channels/generic.jl")
include("channels/onebody.jl")

# circuit
include("circuit.jl")

# qubit hamiltonian
include("qubitsham/qterm.jl")
include("qubitsham/qoperator.jl")

end