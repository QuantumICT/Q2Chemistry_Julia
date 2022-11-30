
"""
	struct QuantumMap{N, T<:Number} <: Gate{N}
"""
struct QuantumMap{N, T<:Number} <: AbstractQuantumMap{N}
	pos::NTuple{N, Int}
	kraus::Vector{Matrix{T}}
	ordered_pos::NTuple{N, Int}
	ordered_data::Matrix{T}

end
function QuantumMap(pos::NTuple{N, Int}, mr::Vector{<:AbstractMatrix{T}}) where {N, T}
	shapes = ntuple(i->2, 2*N)
	m = [reshape(mj, shapes) for mj in mr]
	ordered_pos, ordered_data = _get_norm_order_list(pos, m)
	L = 2^N
	kraus = [convert(Matrix{T}, reshape(mj, L, L)) for mj in m]
	ordered_data = convert(Matrix{T}, _compute_supermat(kraus)) 
	QuantumMap(pos, kraus, ordered_pos, ordered_data)
end


QuantumMap(pos::Vector{Int}, m::Vector{<:AbstractMatrix}) = QuantumMap(Tuple(pos), m)
QuantumMap(pos::Int, m::Vector{<:AbstractMatrix}) = QuantumMap((pos,), m)

positions(x::QuantumMap) = x.pos
kraus_matrices(x::QuantumMap) = x.kraus
ordered_positions(x::QuantumMap) = x.ordered_pos
ordered_supermat(x::QuantumMap) = x.ordered_data
change_positions(x::QuantumMap{N}, m::AbstractDict) where N = QuantumMap(ntuple(i->m[x.pos[i]], N), kraus_matrices(x))



function _get_norm_order_list(key::NTuple{N, Int}, p::Vector{<:AbstractArray}) where N
	seq = sortperm([key...])
	perm = (seq..., [s + N for s in seq]...)
	return key[seq], [permute(item, perm) for item in p]
end


function _compute_supermat(p::Vector{<:AbstractMatrix})
	isempty(p) && error("no kraus operators.")
	d = size(p[1], 1)
	@tensor m[1,3,2,4] := p[1][1,2] * conj(p[1][3,4])
	for i in 2:length(p)
		@tensor m[1,3,2,4] += p[i][1,2] * conj(p[i][3,4])
	end
	return reshape(m, d^2, d^2)
end

