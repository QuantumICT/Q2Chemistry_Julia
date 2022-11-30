
struct QubitsTerm <: QuantumOperation
	positions::Vector{Int}
	op::Vector{<:AbstractMatrix}
	coeff::Number

function QubitsTerm(pos::Vector{Int}, m::Vector, v::Number)
	(length(pos) == length(m)) || error("number of sites mismatch with number of ops.")
	pos, m = _get_normal_order(pos, m)
	new(pos, m, v)
end 

end

positions(x::QubitsTerm) = x.positions
oplist(x::QubitsTerm) = x.op
coeff(x::QubitsTerm) = x.coeff


QubitsTerm(pos::Tuple, m::Vector, v::Number) = QubitsTerm([pos...], m, v)

function QubitsTerm(x::AbstractDict{Int}; coeff::Number=1.)
	sites, ops = dict_to_site_ops(x)
	return QubitsTerm(sites, ops, coeff)
end

QubitsTerm(x::Pair{Int, <:Union{AbstractString, AbstractMatrix}}...; coeff::Number=1.) = QubitsTerm(
	Dict(x...), coeff=coeff)


Base.copy(x::QubitsTerm) = QubitsTerm(copy(positions(x)), copy(oplist(x)), coeff(x))
Base.isempty(x::QubitsTerm) = isempty(oplist(x))
Base.adjoint(x::QubitsTerm) = QubitsTerm(positions(x), [item' for item in oplist(x)], conj(coeff(x)))
Base.:*(x::QubitsTerm, y::Number) = QubitsTerm(positions(x), oplist(x), coeff(x)*y)
Base.:*(x::Number, y::QubitsTerm) = y * x

function Base.eltype(x::QubitsTerm)
	T = typeof(coeff(x))
	for item in oplist(x)
		T = promote_type(T, eltype(item))
	end
	return T
end

function _get_normal_order(key::Vector{Int}, op)
	seq = sortperm(key)
	return key[seq], [_get_op(item) for item in op[seq]]
end


_get_op(m::AbstractMatrix) = m
_get_op(m::String) = _op_mapping[m]

const _op_mapping = Dict("X"=>X, "Y"=>Y, "Z"=>Z, "+"=>[0. 1.; 0. 0.], "-"=>[0. 0.; 1. 0.], 
	"I"=>I₂, "u"=>[1. 0.; 0. 0.], "d"=>[0. 0.; 0. 1.], "0"=>[1. 0.; 0. 0.], "1"=>[0. 0.; 0. 1.])


function dict_to_site_ops(opstr::AbstractDict)
	sites = []
	ops = []
	for (k, v) in opstr
	    push!(sites, k)
	    push!(ops, v)
	end
	return [sites...], [ops...]
end

