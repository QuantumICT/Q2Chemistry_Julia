abstract type AbstractQuantumMap{N} <: QuantumOperation end


# interfaces for Quantum Channel

nqubits(x::AbstractQuantumMap{N}) where N = N

positions(x::AbstractQuantumMap) = x.pos

"""
	kraus_matrices(x::AbstractQuantumMap)
"""
kraus_matrices(x::AbstractQuantumMap) = error("kraus_matrices not implemented for map type $(typeof(x)).")


"""
	ordered_supermat(x::AbstractQuantumMap)
"""
ordered_supermat(x::AbstractQuantumMap) = error("ordered_supermat not implemented for map type $(typeof(x)).")

change_positions(x::AbstractQuantumMap, m::AbstractDict) = error("change_positions not implemented for map type $(typeof(x)).")

shift(x::AbstractQuantumMap, j::Int) = _shift_gate_util!(x, j, Dict{Int, Int}())

is_tp(x::AbstractQuantumMap) = _is_tp(kraus_matrices(x))


function _is_tp(m::Vector{<:AbstractMatrix})
	isempty(m) && error("no kraus operators.")
	out = m[1]' * m[1]
	for i in 2:length(m)
		out .+= m[i]' * m[i]
	end
	return out ≈ I
end

