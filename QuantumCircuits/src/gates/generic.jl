
"""
	struct QuantumGate{N, T<:Number} <: Gate{N}
"""
struct QuantumGate{N, T<:Number} <: Gate{N}
	pos::NTuple{N, Int}
	data::Matrix{T}
	ordered_pos::NTuple{N, Int}
	ordered_data::Matrix{T}


end
function QuantumGate(pos::NTuple{N, Int}, mr::AbstractMatrix{T}) where {N, T}
	m = reshape(mr, ntuple(i->2, 2*N))
	ordered_pos, ordered_data = _get_norm_order(pos, m)
	L = 2^N
	QuantumGate(pos, convert(Matrix{T}, reshape(m, L, L)), ordered_pos, convert(Matrix{T}, reshape(ordered_data, L, L)))
end


QuantumGate(pos::Vector{Int}, m::AbstractMatrix) = QuantumGate(Tuple(pos), m)
QuantumGate(pos::Int, m::AbstractMatrix) = QuantumGate((pos,), m)

positions(x::QuantumGate) = x.pos
mat(x::QuantumGate) = x.data
ordered_positions(x::QuantumGate) = x.ordered_pos
ordered_mat(x::QuantumGate) = x.ordered_data
change_positions(x::QuantumGate{N}, m::AbstractDict) where N = QuantumGate(ntuple(i->m[x.pos[i]], N), mat(x))

gate(pos::Union{NTuple, Vector, Int}, m::AbstractMatrix) = QuantumGate(pos, m)

struct AdjointQuantumGate{N, G<:Gate{N}} <: Gate{N}
	parent::G
end

positions(x::AdjointQuantumGate) = positions(x.parent)
mat(x::AdjointQuantumGate) = mat(x.parent)'
ordered_positions(x::AdjointQuantumGate) = ordered_positions(x.parent)
ordered_mat(x::AdjointQuantumGate) = ordered_mat(x.parent)'
change_positions(x::AdjointQuantumGate, m::AbstractDict) = change_positions(x.parent, m)'


Base.adjoint(x::Gate) = AdjointQuantumGate(x)
Base.adjoint(x::AdjointQuantumGate) = x.parent


