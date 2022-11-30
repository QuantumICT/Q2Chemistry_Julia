

_is_commutative(a::Gate, b::Gate) = _is_commutative(ordered_positions(a), ordered_positions(b))
function _is_commutative(a::Vector{<:Gate}, b::Gate)
    for item in a
        _is_commutative(item, b) || return false
    end
    return true
end
_is_commutative(a::QCircuit, b::Gate) = _is_commutative(a.operations, b)

function commutative_blocks(s::QCircuit)
    t = Vector{typeof(s)}()
    length(s)==0 && return t
    tmp = similar(s)
    for gate in s
        if _is_commutative(tmp, gate)
            push!(tmp, gate)
        else
            push!(t, tmp)
            tmp = similar(s)
            push!(tmp, gate)
        end
    end
    isempty(tmp) && error("something wrong.")
    push!(t, tmp)
    return t
end


function try_mul_two_ops(a::Gate, b::Gate, workspace::AbstractVector)
    a_op = ordered_op(a)
    b_op = ordered_op(b)
    L = length(a_op) + length(b_op)
    if length(workspace) < L
        resize!(workspace, L)
    end
    if !isa(a_op, StridedArray)
        a_op = copyto!(reshape(view(workspace, 1:length(a_op)), size(a_op)), a_op)
    end
    if !isa(b_op, StridedArray)
       b_op = copyto!(reshape(view(workspace, length(a_op)+1:L), size(b_op)), b_op) 
    end
    r = QuantumSpins._try_mul_two_ops_impl(ordered_positions(a), a_op, ordered_positions(b), b_op)
    (r === nothing) && return nothing
    k = (length(QuantumCircuits.positions(a)) >= length(QuantumCircuits.positions(b))) ? QuantumCircuits.positions(a) : QuantumCircuits.positions(b)
    _lk = 2^(length(k))
    return gate(k, reshape(r, _lk, _lk))
end

function try_absorb_one_gate_by_mul!(a::Vector, gt::Gate, workspace::AbstractVector)
    for i in length(a):-1:1
        r = try_mul_two_ops(gt, a[i], workspace)
        if r !== nothing
            a[i] = r
            return true
        end
        _is_commutative(a[i], gt) || return false
    end
    return false
end

function try_absorb_by_mul_impl(a::Vector, workspace::AbstractVector)
    isempty(a) && return a
    L = length(a)
    b = copy(a)
    for i in L:-1:2
        item = b[i]
        br = b[1:(i-1)]
        if try_absorb_one_gate_by_mul!(br, item, workspace)
            b = [br; b[(i+1):end]]
        end
    end
    return b
end

_scalar_t(a::Gate) = eltype(QuantumCircuits.mat(a))

function _scalar_t(a::QCircuit)
	T = Float64
	for g in a
		T = promote_type(T, _scalar_t(g) )
	end
	return T
end

function _fuse_gate_impl(a::QCircuit)
    workspace = _scalar_t(a)[]
    b = try_absorb_by_mul_impl(a.operations, workspace)
    r = similar(a)
    append!(r, b)
    return r
end

function _flatten(a::QCircuit)
	r = similar(a)
	for item in a
		if isa(item, Gate)
			push!(r, item)
		elseif isa(item, QCircuit)
			append!(r, _flatten(item))
		else
			error("unexpected type $(typeof(item))")
		end
	end
	return r
end

function _fuse_gates(a::QCircuit)
    b = _fuse_gate_impl(a')
    return _fuse_gate_impl(b')
end
fuse_gates(a::QCircuit) = _fuse_gates(_flatten(a))

