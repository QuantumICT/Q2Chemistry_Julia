

# one body gates
struct XGate <: Gate{1}
	pos::Tuple{Int}
end
XGate(pos::Int) = XGate((pos,))

mat(s::XGate) = X
change_positions(s::XGate, m::AbstractDict) = XGate((m[s.pos[1]],))

Base.adjoint(x::XGate) = x


struct YGate <: Gate{1}
	pos::Tuple{Int}
end
YGate(pos::Int) = YGate((pos,))

mat(x::YGate) = Y
change_positions(s::YGate, m::AbstractDict) = YGate((m[s.pos[1]],))


struct ZGate <: Gate{1}
	pos::Tuple{Int}
end
ZGate(pos::Int) = ZGate((pos,))

mat(x::ZGate) = Z
change_positions(s::ZGate, m::AbstractDict) = ZGate((m[s.pos[1]],))


struct SGate <: Gate{1}
	pos::Tuple{Int}
end
SGate(pos::Int) = SGate((pos,))

mat(x::SGate) = S
change_positions(s::SGate, m::AbstractDict) = SGate((m[s.pos[1]],))


struct sqrtXGate <: Gate{1}
	pos::Tuple{Int}
end
sqrtXGate(pos::Int) = sqrtXGate((pos,))

mat(x::sqrtXGate) = sqrtX
change_positions(s::sqrtXGate, m::AbstractDict) = sqrtXGate((m[s.pos[1]],))


struct sqrtYGate <: Gate{1}
	pos::Tuple{Int}
end
sqrtYGate(pos::Int) = sqrtYGate((pos,))

mat(x::sqrtYGate) = sqrtY
change_positions(s::sqrtYGate, m::AbstractDict) = sqrtYGate((m[s.pos[1]],))


struct HGate <: Gate{1}
	pos::Tuple{Int}
end
HGate(pos::Int) = HGate((pos,))

mat(x::HGate) = H
change_positions(s::HGate, m::AbstractDict) = HGate((m[s.pos[1]],))

struct TGate <: Gate{1}
	pos::Tuple{Int}
end
TGate(pos::Int) = TGate((pos,))

mat(x::TGate) = T
change_positions(s::TGate, m::AbstractDict) = TGate((m[s.pos[1]],))


# parameteric one body gate
struct RxGate <: ParametricGate{1}
	pos::Tuple{Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
RxGate(pos::Int, p::Real; isparas::Bool=false) = RxGate((pos,), [convert(Float64, p)], [isparas])
RxGate(pos::Int; θ::Real, isparas::Bool=false) = RxGate(pos, θ, isparas=isparas)

mat(x::RxGate) = Rx(x.paras[1])
change_positions(s::RxGate, m::AbstractDict) = RxGate((m[s.pos[1]],), s.paras, s.isparas)


struct RyGate <: ParametricGate{1}
	pos::Tuple{Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
RyGate(pos::Int, p::Real; isparas::Bool=false) = RyGate((pos,), [convert(Float64, p)], [isparas])
RyGate(pos::Int; θ::Real, isparas::Bool=false) = RyGate(pos, θ, isparas=isparas)

mat(x::RyGate) = Ry(x.paras[1])
change_positions(s::RyGate, m::AbstractDict) = RyGate((m[s.pos[1]],), s.paras, s.isparas)


struct RzGate <: ParametricGate{1}
	pos::Tuple{Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
RzGate(pos::Int, p::Real; isparas::Bool=false) = RzGate((pos,), [convert(Float64, p)], [isparas])
RzGate(pos::Int; θ::Real, isparas::Bool=false) = RzGate(pos, θ, isparas=isparas)

mat(x::RzGate) = Rz(x.paras[1])
change_positions(s::RzGate, m::AbstractDict) = RzGate((m[s.pos[1]],), s.paras, s.isparas)


struct PHASEGate <: ParametricGate{1}
	pos::Tuple{Int}
	paras::Vector{Float64}
	isparas::Vector{Bool}
end
PHASEGate(pos::Int, p::Real; isparas::Bool=false) = PHASEGate((pos,), [convert(Float64, p)], [isparas])
PHASEGate(pos::Int; θ::Real, isparas::Bool=false) = PHASEGate(pos, θ, isparas=isparas)

mat(x::PHASEGate) = PHASE(x.paras[1])
change_positions(s::PHASEGate, m::AbstractDict) = PHASEGate((m[s.pos[1]],), s.paras, s.isparas)

Base.adjoint(x::PHASEGate) = PHASEGate(positions(x), -x.paras, x.isparas)

