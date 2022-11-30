

struct TOFFOLIGate <: Gate{3}
	pos::Tuple{Int, Int, Int}
	ordered_pos::Tuple{Int, Int, Int}
	ordered_data::Matrix{Float64}
end

function TOFFOLIGate(pos::Tuple{Int, Int, Int})
	m = reshape(TOFFOLI, ntuple(i->2, 6))
	ordered_pos, ordered_data = _get_norm_order(pos, m)
	return TOFFOLIGate(pos, ordered_pos, reshape(ordered_data, 8, 8))
end
TOFFOLIGate(a::Int, b::Int, c::Int) = TOFFOLIGate((a,b,c))

mat(x::TOFFOLIGate) = TOFFOLI
ordered_positions(x::TOFFOLIGate) = x.ordered_pos
ordered_mat(x::TOFFOLIGate) = x.ordered_data
change_positions(x::TOFFOLIGate, m::AbstractDict) = TOFFOLIGate(ntuple(i->m[x.pos[i]], 3))

Base.adjoint(x::TOFFOLIGate) = x




struct FREDKINGate <: Gate{3}
	pos::Tuple{Int, Int, Int}
	ordered_pos::Tuple{Int, Int, Int}
	ordered_data::Matrix{Float64}
end

function FREDKINGate(pos::Tuple{Int, Int, Int})
	m = reshape(FREDKIN, ntuple(i->2, 6))
	ordered_pos, ordered_data = _get_norm_order(pos, m)
	return FREDKINGate(pos, ordered_pos, reshape(ordered_data, 8, 8))

end
FREDKINGate(a::Int, b::Int, c::Int) = FREDKINGate((a, b, c))


mat(x::FREDKINGate) = FREDKIN
ordered_positions(x::FREDKINGate) = x.ordered_pos
ordered_mat(x::FREDKINGate) = x.ordered_data
change_positions(x::FREDKINGate, m::AbstractDict) = FREDKINGate(ntuple(i->m[x.pos[i]], 3))

Base.adjoint(x::FREDKINGate) = x


# control control gate
struct CONTROLCONTROLGate{M<:AbstractMatrix} <: Gate{3}
	pos::Tuple{Int, Int, Int}
	target::M
end
CONTROLCONTROLGate(a::Int, b::Int, c::Int, m::AbstractMatrix) = CONTROLCONTROLGate((a,b,c), m) 
mat(x::CONTROLCONTROLGate) = CONTROLCONTROL(x.target)

change_positions(x::CONTROLCONTROLGate, m::AbstractDict) = CONTROLCONTROLGate(ntuple(i->m[x.pos[i]], 3), x.target)


struct CCPHASEGate <: ParametricGate{3}
	pos::Tuple{Int, Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
CCPHASEGate(key::Tuple{Int, Int, Int}, p::Real; isparas::Bool=true) = CCPHASEGate(key, [convert(Float64, p)], [isparas])
# CCPHASEGate(key::Tuple{Int, Int, Int}; θ::Real, isparas::Bool=true) = CCPHASEGate(key, θ, isparas=isparas)
CCPHASEGate(a::Int, b::Int, c::Int, θ::Real, isparas::Bool=true) = CCPHASEGate((a, b, c), θ, isparas=isparas)


mat(x::CCPHASEGate) = CONTROLCONTROL(PHASE(x.paras[1]))
change_positions(x::CCPHASEGate, m::AbstractDict) = CCPHASEGate(ntuple(i->m[x.pos[i]], 3), x.paras, x.isparas)
Base.adjoint(x::CCPHASEGate) = CCPHASEGate(positions(x), -x.paras, x.isparas)

