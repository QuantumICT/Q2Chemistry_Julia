abstract type QuantumCircuit <: QuantumOperation end

"""
	struct QMeasure <: QuantumOperation
	Quantum Measurement
"""
struct QMeasure <: QuantumOperation
	position::Int
	auto_reset::Bool
	keep::Bool

	QMeasure(key::Int; auto_reset::Bool=true, keep::Bool=true) = new(key, auto_reset, keep)
end


"""
	struct QSelect <: QuantumOperation
	post selection, do we need it?
"""
struct QSelect <: QuantumOperation
	position::Int
	state::Int

	function QSelect(key::Int, state::Int=0)
		(state==0 || state==1) || error("selected state must be 0 or 1.")
		new(key, state)
	end
end


const CircuitElement = Union{QuantumCircuit, Gate, AbstractQuantumMap}


"""
	struct QCircuit <: QuantumOperation
	Quantum Circuit
"""
struct QCircuit <: QuantumCircuit
	operations::Vector{CircuitElement}
end
QCircuit() = QCircuit(Vector{CircuitElement}())
QCircuit(x::QCircuit) = QCircuit(copy(x.operations))
Base.similar(circuit::QCircuit) = QCircuit()
Base.copy(circuit::QCircuit) = QCircuit(copy(circuit.operations))
Base.length(x::QCircuit) = length(x.operations)
Base.getindex(x::QCircuit, j::Int) = x.operations[j]
Base.setindex!(x::QCircuit, v, j::Int) = setindex!(x.operations, v, j)
Base.iterate(x::QCircuit) = iterate(x.operations)
Base.iterate(x::QCircuit, state) = iterate(x.operations, state)
Base.eltype(x::QCircuit) = eltype(x.operations)
Base.empty!(x::QCircuit) = empty!(x.operations)
Base.isempty(x::QCircuit) = isempty(x.operations)
Base.firstindex(x::QCircuit) = firstindex(x.operations)
Base.lastindex(x::QCircuit) = lastindex(x.operations)
Base.reverse(x::QCircuit) = QCircuit(reverse(x.operations))

"""
	Base.:*(x::QCircuit, y::QCircuit)
	Do we need it?
"""
Base.:*(x::QCircuit, y::QCircuit) = QCircuit(vcat(y.operations, x.operations))
Base.adjoint(x::QCircuit) = QCircuit(adjoint.(Iterators.reverse(x.operations)))

Base.push!(x::QCircuit, s::Gate) = push!(x.operations, s)
Base.push!(x::QCircuit, s::QCircuit) = push!(x.operations, s)
Base.append!(x::QCircuit, s::Vector{<:CircuitElement}) = append!(x.operations, s)
Base.append!(x::QCircuit, s::QCircuit) = append!(x.operations, s.operations)
Base.push!(x::QCircuit, s::Tuple{Int, <:AbstractMatrix}) = push!(x, gate(s[1], s[2]))
Base.push!(x::QCircuit, s::Tuple{NTuple{N, Int}, <:AbstractMatrix}) where N = push!(x, gate(s[1], s[2]))

change_positions(x::QCircuit, m::AbstractDict) = QCircuit([change_positions(o, m) for o in x.operations])

shift(x::QCircuit, j::Int) = _shift_gate_util!(x, j, Dict{Int, Int}())

function positions(x::QCircuit)
	pos = Set{Int}()
	for g in x
		r = positions(g)
		for item in r
			push!(pos, item)
		end
	end
	return sort([pos...])
end

function nparameters(x::QCircuit)
	n = 0
	for o in x
		n += nparameters(o)
	end
	return n
end

function parameters(x::QCircuit)
	paras = Float64[]
	for o in x
		p = parameters(o)
		if !isnothing(p)
			append!(paras, p)
		end
	end
	return paras
end
function active_parameters(x::QCircuit)
	paras = Float64[]
	for o in x
		p = active_parameters(o)
		if !isnothing(p)
			append!(paras, p)
		end
	end
	return paras
end
function activate_parameter!(x::QCircuit, j::Int)
	activate_parameters!(x[j])
	return x
end
function activate_parameters!(x::QCircuit)
	for o in x
		activate_parameters!(o)
	end
	return x
end
function deactivate_parameter!(x::QCircuit, j::Int)
	deactivate_parameters!(x[j])
	return x
end
function deactivate_parameters!(x::QCircuit)
	for o in x
		deactivate_parameters!(o)
	end
	return x
end

function reset_parameters!(x::QCircuit, p::Vector{<:Real})
	(nparameters(x) == length(p)) || throw("number of parameters mismatch.")
	pos = reset_parameters_util!(x, p, 0)
	@assert pos == length(p)
	return x
end


function reset_parameters_util!(x::QCircuit, p::Vector{<:Real}, pos::Int)
	for o in x
		pos = reset_parameters_util!(o, p, pos)
	end
	return pos
end

function _shift_gate_util!(x::QCircuit, j::Int, m::AbstractDict)
	r = similar(x)
	for item in x
		push!(r, _shift_gate_util!(item, j, m))
	end
	return r
end
