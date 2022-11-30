using MPI
MPI.Init()
const root = 0
using LinearAlgebra
push!(LOAD_PATH, "./QuantumSpins/src")
push!(LOAD_PATH, "./QuantumCircuits/src")
push!(LOAD_PATH, "./MPSSimulator/src")
using QuantumSpins
using QuantumCircuits
using MPSSimulator
using NLopt


if length(ARGS)>1
    D_global = parse(Int, ARGS[2])
else
    D_global = 128
end
if MPI.Comm_rank(MPI.COMM_WORLD) == root
    println("Bond dimension = $D_global")
end

const QubitOperator = Vector{Pair{Tuple{Vararg{Tuple{Int64, String}}}, ComplexF64}}
const Hamiltonian = QuantumCircuits.QubitsOperator
function Hamiltonian(x::QubitOperator)
	qop = Hamiltonian()
	for (k, v) in x
		d = Dict{Int, Union{AbstractString, AbstractMatrix}}()
		for (idx, qubit_gate) in k
			d[idx] = qubit_gate
		end
		if d.count == 0
			d[1] = "I"
		end
		QuantumCircuits.add!(qop, QubitsTerm(d, coeff=v))
	end
	return qop	
end


const hy_matrix = [+1.0+0.0im +0.0-1.0im; 0.0+1.0im -1.0+0.0im] * 0.5 * 2^0.5
function p_range(n::Int, n_procs::Int, pid::Int)::UnitRange{Int}
	aprocs = n_procs - n % n_procs + 1
	q = n ÷ n_procs
	if (pid < aprocs)
			pstart = (pid - 1) * q + 1
			pend = pstart + q - 1
	else
			pstart = (aprocs-1) * q + (pid - aprocs)*(q+1) + 1
			pend = pstart + q
	end
	return pstart:pend
end

function MPI.Scatterv(QO :: Union{Nothing, Hamiltonian}, root :: Integer, comm::MPI.Comm) :: Hamiltonian
	if root == MPI.Comm_rank(comm)
			@assert QO != nothing
			len = QO.data.count
			n_procs = MPI.Comm_size(comm)
			QO_pairs = collect(QO.data)
			reqs = Vector{MPI.Request}(undef, n_procs - 1)
			for i in 1:n_procs-1
					qo = Hamiltonian(QuantumCircuits.QOP_DATA_TYPE(@view QO_pairs[p_range(len, n_procs, i+1)]))
					reqs[i] = MPI.isend(qo, i, 0, comm)
			end
			qo_self = Hamiltonian(QuantumCircuits.QOP_DATA_TYPE(@view QO_pairs[p_range(len, n_procs, 1)]))
			
			MPI.Waitall!(reqs)
			return qo_self
	else
			qo, status = MPI.recv(root, 0, comm)
			return qo
	end
end

# Return circuit and map
function construct_vqe_circuit(occ_indices_spin :: Vector{Int}, ucc_operator_pool_qubit_op::Vector)
	circuit = QCircuit()
	m = Vector{Tuple{Float64, Int}}()
	for i in occ_indices_spin
		push!(circuit, XGate(i))
	end
	for i in 1:length(ucc_operator_pool_qubit_op)

		for (terms,coeff) in ucc_operator_pool_qubit_op[i]
			length(terms) == 0 && continue
			for (idx, qubit_gate) in terms
				if qubit_gate == "X"
					push!(circuit, HGate(idx))
				elseif qubit_gate == "Y"
					push!(circuit, gate(idx, hy_matrix))
				end
			end

			for i in length(terms)-1:-1:1
				push!(circuit, CNOTGate(terms[i+1][1], terms[i][1]))
			end

			push!(circuit, RzGate(terms[1][1], 0.0, isparas=true))
			push!(m, (-2.0*imag(coeff),i))

			for i in 1:length(terms)-1
				push!(circuit, CNOTGate(terms[i+1][1], terms[i][1]))
			end

			for (idx, qubit_gate) in terms
				if qubit_gate == "X"
					push!(circuit, HGate(idx))
				elseif qubit_gate == "Y"
					push!(circuit, gate(idx, hy_matrix))
				end
			end

		end

	end
	return circuit, m
end

function reset_parameters_wrapper(p::Vector{<:Real}; x::QCircuit) :: QCircuit
	reset_parameters!(x, p)
	return x
end

function apply_map(amps :: Vector{Float64}; map :: Vector{Tuple{Float64, Int}}) :: Vector{Float64}
	return [line[1]*amps[line[2]] for line in map]
end

# Normalization of MPS
normalize_mps(psi::MPS) = normalize(psi)


function vqe(parameters :: Dict, comm :: MPI.Comm)

	rank = MPI.Comm_rank(comm)
	n_qubits::Int = parameters["n_qubits"]
	n_amplitudes::Int = parameters["n_params"]
	observable::Hamiltonian = parameters["hamiltonian_qubit_op"]
	spin_orbital_occupied_indices::Vector{Int} = parameters["spin_orbital_occupied_indices"]
	ucc_operator_pool_qubit_op::Vector{QubitOperator} = parameters["ucc_operator_pool_qubit_op"]	
	circuit, map = construct_vqe_circuit(spin_orbital_occupied_indices, ucc_operator_pool_qubit_op)
	initial_state = statevector_mps(n_qubits)
	trunc = MPSTruncation(ϵ=1.0e-6, D=D_global)
	loss(amplitudes) = begin
		temp = normalize_mps(apply(reset_parameters_wrapper(apply_map(amplitudes, map=map), x=circuit), initial_state, trunc=trunc))
		real(expectation(observable, temp, trunc=trunc))
	end
	

	if rank == root
		println("qubits: ", n_qubits)
		println("parameters: ",length(ucc_operator_pool_qubit_op))
		println("parametric Rz gates: ",length(map))
		println("local hamiltonian terms: ",observable.data.count)


		function f(x::Vector, grad::Vector)
            MPI.bcast(x, root, comm)
                v_loc = loss(x)
                v = MPI.Reduce(v_loc, +, root, comm)
	            println("Energy: $v")
            v
		end
		work = () -> begin

                opt = Opt(:LN_BOBYQA, n_amplitudes)
                opt.min_objective = f
                opt.maxeval = 1000
                opt.ftol_abs = 1e-8

            (minf, minx, ret) = optimize(opt, fill(0.0, n_amplitudes))

			MPI.bcast(nothing, root, comm)
			println("Optimized energy result: $(minf)")
		end

	else
        
		    work = () -> begin
			    while (true)
				    x0::Union{Nothing, Vector{Float64}} = MPI.bcast(nothing, root, comm)
				    if x0 == nothing
					    break
				    end
                    e_loc = loss(x0)
				    MPI.Reduce(e_loc, +, root, comm)
			    end		
		    end

	end


	work()


end

function read_binary_qubit_op(f::IO, n_qubits::Int) :: QubitOperator

    pauli_symbol_dict = Dict(
        0 => "I",
        1 => "X",
        2 => "Y",
        3 => "Z"
    )
	qubit_op_dict = QubitOperator()
	len = read(f, Int64)
	if len > 0			# paulis are stored in dense array
		for i in 1:len
			coeff_tmp = read(f, ComplexF64)
			pauli_str_tmp = Int8[read(f, Int8) for i in 1:n_qubits]
			pauli_str_tuple = Tuple([(i, pauli_symbol_dict[pauli_str_tmp[i]])
                                 for i in 1:n_qubits
                                 if pauli_str_tmp[i] != 0])
			push!(qubit_op_dict, pauli_str_tuple => coeff_tmp)
		end
	else				# paulis are stored in compressed array
		len = -len
		for i in 1:len
			coeff_tmp = read(f, ComplexF64)
			len_pauli_array = read(f, Int32)
			pauli_str_tuple = Vector{Tuple{Int64, String}}()
			for i in 1:len_pauli_array
                # Indices are shifted by 1 due to array indices start at 1 in Julia.
				push!(pauli_str_tuple, (read(f, Int32) + 1, pauli_symbol_dict[read(f, Int8)]))
			end
			push!(qubit_op_dict, Tuple(pauli_str_tuple) => coeff_tmp)
		end
	end
	return qubit_op_dict
end

function read_binary_dict(f::IO)::Dict{String,Any}

    identifier = read(f, Float64)
    err_msg = identifier != Float64(99.0212) && error("The file is not saved as vqe parameters.")

	d = Dict{String, Any}()
	d["n_qubits"] = Int64(read(f, Int32))
	len_spin_indices = read(f, Int32)

    # Indices are shifted by 1 due to array indices start at 1 in Julia.
	d["spin_orbital_occupied_indices"] = Int64[read(f,Int32) + 1 for i in 1:len_spin_indices]

	d["n_params"] = Int64(read(f, Int32))
	d["ucc_operator_pool_qubit_op"] = Vector{QubitOperator}(undef, d["n_params"])
	for i in 1:d["n_params"]
		d["ucc_operator_pool_qubit_op"][i] = read_binary_qubit_op(f, d["n_qubits"])
	end
	return d
end
xzqsize(x) = Base.format_bytes(Base.summarysize(x))

function main()

	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	p = Dict{String, Any}()
	if rank == root
		f = open(ARGS[1], "r")
		p = read_binary_dict(f)
		MPI.bcast(p, root, comm)


		q = read_binary_qubit_op(f, p["n_qubits"])
		println("original global hamiltonian terms: ",length(q))


		H = Hamiltonian(q)
		println("global hamiltonian terms: ",H.data.count)

		
		p["hamiltonian_qubit_op"] = MPI.Scatterv(H, root, comm)
	else
		p = MPI.bcast(nothing, root, comm)
		p["hamiltonian_qubit_op"] = MPI.Scatterv(nothing, root, comm)
	end

	vqe(p, comm)
end
main()
