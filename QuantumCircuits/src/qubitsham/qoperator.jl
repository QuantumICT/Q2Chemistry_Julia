const QOP_DATA_TYPE = Dict{Tuple{Int, Vararg{Int, N} where N}, Vector{Tuple{Vector{AbstractMatrix}, Number}}}
const QOP_DATA_VALUE_TYPE = Vector{Tuple{Vector{AbstractMatrix}, Number}}

struct QubitsOperator <: QuantumOperation
	data::QOP_DATA_TYPE
end

QubitsOperator() = QubitsOperator(QOP_DATA_TYPE())

Base.copy(x::QubitsOperator) = QubitsOperator(copy(x.data))
Base.keys(x::QubitsOperator) = keys(x.data)

function QubitsOperator(x::QubitsTerm, v::QubitsTerm...)
	r = QubitsOperator()
	add!(r, x)
	for item in v
	    add!(r, item)
	end
	return r
end 

function absorb_one_bodies(ham::QubitsOperator)
	r = QubitsOperator()
	# L = get_largest_pos(ham)
	iden = I₂
	for (key, value) in ham.data
		if length(key)==1
			i = key[1]
			for (k2, v2) in ham.data
			    if ((length(k2) > 1) && (i in k2))
			    	for item in value
			    		ops = [ (kk == i) ? item[1][1] : iden for kk in k2]
			    		add!(r, QubitsTerm([k2...], ops, item[2]))
			        end
			        break
			    end
			end
		else
			for item in value
			    add!(r, QubitsTerm([key...], item[1], item[2]))
			end
		end
	end
	return r
end

simplify(h::QubitsOperator) = absorb_one_bodies(h)

function add!(x::QubitsOperator, m::QubitsTerm)
	pos = Tuple(positions(m))
	v = get!(x.data, pos, QOP_DATA_VALUE_TYPE())
	push!(v, convert(Tuple{Vector{AbstractMatrix}, Number}, (oplist(m), coeff(m))) )
	return x
end

function Base.:+(x::QubitsTerm, y::QubitsTerm)
	r = QubitsOperator()
	add!(r, x)
	add!(r, y)
	return r
end
function Base.:+(x::QubitsTerm, y::QubitsOperator)
	r = copy(y)
	add!(r, x)
	return r
end
Base.:+(x::QubitsOperator, y::QubitsTerm) = y + x
function Base.:+(x::QubitsOperator, y::QubitsOperator)
	r = copy(x)
	for (k, v) in y.data
		kk = [k...]
		for item in v
		    add!(r, QubitsTerm(kk, item[1], item[2]))
		end
	end
	return r
end


function matrix(L::Int, x::QubitsOperator)
	# is_constant(x) || error("input must be a constant operator.")
	h = nothing
	for (k, v) in x.data
		for item in v
			tmp = _generateprodham(L, sites_ops_to_dict(k, item[1])) * item[2]
			if isnothing(h)
			    h = tmp
			else
				h += tmp
			end
		end
	end
	return h
end

matrix(x::QubitsOperator) = matrix(get_largest_pos(x), x)
matrix(L::Int, x::QubitsTerm) = _generateprodham(L, sites_ops_to_dict(positions(x), oplist(x))) * coeff(x)


change_positions(x::QubitsOperator, m::AbstractDict) = QubitsOperator(QOP_DATA_TYPE(_index_map(k, m)=>v for (k, v) in x.data)) 

function Base.eltype(x::QubitsOperator)
	T = Int
	for (k, v) in x.data
		for (o, c) in v
			T = promote_type(typeof(c), T)
			for oj in o
				T = promote_type(T, eltype(oj))
			end
		end
	end
	return T
end

Base.adjoint(x::QubitsOperator) = QubitsOperator(QOP_DATA_TYPE(k=>[(adjoint.(a), conj(b)) for (a, b) in v] for (k, v) in x.data))
Base.:*(x::QubitsOperator, y::Number) = QubitsOperator(QOP_DATA_TYPE(k=>[(a, b*y) for (a, b) in v] for (k, v) in x.data))
Base.:*(y::Number, x::QubitsOperator) = x * y

# function Base.getindex(h::QubitsOperator, v::Vector{Int})
# 	index_map = Dict(vj=>j for (j, vj) in enumerate(v))
# 	data_new = typeof(data(h))()
# 	for (k, bond) in data(h)
# 	    if issubset(k, v)
# 	    	new_k = Tuple(index_map[item] for item in k)
# 	    	seq = sortperm([new_k...])
# 	    	new_k = new_k[seq]
# 	    	data_new[new_k] = permute(bond, seq)
# 	    end
# 	end
# 	return QubitsOperator(data_new)
# end

function _generateprodham(L::Int, opstr::Dict{Int, <:AbstractMatrix}) 
	(max(keys(opstr)...) <= L) || error("op str out of bounds")
	i = 1
	ops = []
	for i in 1:L
	    v = get(opstr, i, nothing)
	    if v === nothing
	        v = I₂
	    else
	    	(size(v, 1)==2 && size(v, 2)==2) || error("dimension mismatch with dim.")
	    end
	    push!(ops, sparse(v))
	end
	return _kron_ops(reverse(ops))
end

_index_map(x::NTuple{N, Int}, mm::AbstractDict) where N = ntuple(j -> mm[x[j]], N)



function get_largest_pos(x::QubitsOperator)
	L = 0
	for (k, v) in x.data
	    L = max(L, maximum(k))
	end
	return L
end


sites_ops_to_dict(sites, ops)= Dict(sites[j]=>ops[j] for j=1:length(sites))
