
struct SWAPGate <: Gate{2}
	pos::Tuple{Int, Int}

	function SWAPGate(key::Tuple{Int, Int})
		new(Tuple(sort([key...])))
	end
end
SWAPGate(c::Int, t::Int) = SWAPGate((c, t))

mat(x::SWAPGate) = SWAP
ordered_mat(x::SWAPGate) = SWAP
change_positions(x::SWAPGate, m::AbstractDict) = SWAPGate(ntuple(i->m[x.pos[i]], 2))

Base.adjoint(x::SWAPGate) = x


struct iSWAPGate <: Gate{2}
	pos::Tuple{Int, Int}

	function iSWAPGate(key::Tuple{Int, Int})
		new(Tuple(sort([key...])))
	end
end
iSWAPGate(c::Int, t::Int) = iSWAPGate((c, t))

mat(x::iSWAPGate) = iSWAP
ordered_mat(x::iSWAPGate) = iSWAP
change_positions(x::iSWAPGate, m::AbstractDict) = iSWAPGate(ntuple(i->m[x.pos[i]], 2))



# # two body controled gates
# function _get_norm_order(key::NTuple{N, Int}) where N
# 	seq = sortperm([key...])
# 	perm = [seq; [s + N for s in seq]]
# 	return key[seq], perm
# end


struct CZGate <: Gate{2}
	pos::Tuple{Int, Int}

	function CZGate(key::Tuple{Int, Int})
		new(Tuple(sort([key...])))
	end
end
CZGate(c::Int, t::Int) = CZGate((c, t))
CZGate(;control::Int, target::Int) = CZGate(control, target)

mat(x::CZGate) = CZ
ordered_mat(x::CZGate) = CZ
change_positions(x::CZGate, m::AbstractDict) = CZGate(ntuple(i->m[x.pos[i]], 2))

Base.adjoint(x::CZGate) = x


struct CNOTGate <: Gate{2}
	pos::Tuple{Int, Int}
end
CNOTGate(c::Int, t::Int) = CNOTGate((c, t))
CNOTGate(;control::Int, target::Int) = CNOTGate(control, target)

mat(x::CNOTGate) = CNOT
change_positions(x::CNOTGate, m::AbstractDict) = CNOTGate(ntuple(i->m[x.pos[i]], 2))

Base.adjoint(x::CNOTGate) = x


# control gate
struct CONTROLGate{M<:AbstractMatrix} <: Gate{2}
	pos::Tuple{Int, Int}
	target::M
end
CONTROLGate(a::Int, b::Int, target::AbstractMatrix) = CONTROLGate((a, b), target)
CONTROLGate(m::AbstractMatrix; control::Int, target::Int) = CONTROLGate(control, target, m)
mat(x::CONTROLGate) = CONTROL(x.target)

change_positions(x::CONTROLGate, m::AbstractDict) = CONTROLGate(ntuple(i->m[x.pos[i]], 2), x.target)


# parameteric two body gates
struct CRxGate <:ParametricGate{2}
	pos::Tuple{Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
CRxGate(key::Tuple{Int, Int}, p::Real; isparas::Bool=false) = CRxGate(key, [convert(Float64, p)], [isparas])
CRxGate(c::Int, t::Int, p::Real; isparas::Bool=false) = CRxGate((c, t), p, isparas=isparas)
CRxGate(;control::Int, target::Int, θ::Real, isparas::Bool=false) = CRxGate(control, target, θ, isparas=isparas)

mat(x::CRxGate) = CONTROL(Rx(x.paras[1]))
change_positions(x::CRxGate, m::AbstractDict) = CRxGate(ntuple(i->m[x.pos[i]], 2), x.paras, x.isparas)



struct CRyGate <: ParametricGate{2}
	pos::Tuple{Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
CRyGate(key::Tuple{Int, Int}, p::Real; isparas::Bool=false) = CRyGate(key, [convert(Float64, p)], [isparas])
CRyGate(c::Int, t::Int, p::Real; isparas::Bool=false) = CRyGate((c, t), p, isparas=isparas)
CRyGate(;control::Int, target::Int, θ::Real, isparas::Bool=false) = CRyGate(control, target, θ, isparas=isparas)

mat(x::CRyGate) = CONTROL(Ry(x.paras[1]))
change_positions(x::CRyGate, m::AbstractDict) = CRyGate(ntuple(i->m[x.pos[i]], 2), x.paras, x.isparas)


struct CRzGate <: ParametricGate{2}
	pos::Tuple{Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
CRzGate(key::Tuple{Int, Int}, p::Real; isparas::Bool=false) = CRzGate(key, [convert(Float64, p)], [isparas])
CRzGate(c::Int, t::Int, p::Real; isparas::Bool=false) = CRzGate((c, t), p, isparas=isparas)
CRzGate(;control::Int, target::Int, θ::Real, isparas::Bool=false) = CRzGate(control, target, θ, isparas=isparas)

mat(x::CRzGate) = CONTROL(Rz(x.paras[1]))
change_positions(x::CRzGate, m::AbstractDict) = CRzGate(ntuple(i->m[x.pos[i]], 2), x.paras, x.isparas)


struct CPHASEGate <: ParametricGate{2}
	pos::Tuple{Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
CPHASEGate(key::Tuple{Int, Int}, p::Real; isparas::Bool=false) = CPHASEGate(key, [convert(Float64, p)], [isparas])
CPHASEGate(c::Int, t::Int, p::Real; isparas::Bool=false) = CPHASEGate((c, t), p, isparas=isparas)
CPHASEGate(;control::Int, target::Int, θ::Real, isparas::Bool=false) = CPHASEGate(control, target, θ, isparas=isparas)

mat(x::CPHASEGate) = CONTROL(PHASE(x.paras[1]))
change_positions(x::CPHASEGate, m::AbstractDict) = CPHASEGate(ntuple(i->m[x.pos[i]], 2), x.paras, x.isparas)
Base.adjoint(x::CPHASEGate) = CPHASEGate(positions(x), -x.paras, x.isparas)


# theta, phi, deltap, deltam, deltamoff
struct FSIMGate <: ParametricGate{2}
	pos::Tuple{Int, Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
function FSIMGate(key::Tuple{Int, Int}, p::Vector{<:Real}; isparas::Vector{Bool}=[false for j in p])
	(length(p) == 5) || throw(ArgumentError("5 parameters expected."))
	return FSIMGate(key, convert(Vector{Float64}, p) , isparas)
end 
FSIMGate(c::Int, t::Int, p::Vector{<:Real}; isparas::Vector{Bool}=[false for j in p]) = FSIMGate((c, t), p, isparas=isparas)
FSIMGate(;control::Int, target::Int, θs::Vector{<:Real}, isparas::Vector{Bool}=[false for j in θs]) = FSIMGate(control, target, θs, isparas=isparas)

mat(x::FSIMGate) = FSIM(x.paras...)
change_positions(x::FSIMGate, m::AbstractDict) = FSIMGate(ntuple(i->m[x.pos[i]], 2), x.paras, x.isparas)


