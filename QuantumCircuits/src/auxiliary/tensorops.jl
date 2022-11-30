
rkron(a::AbstractMatrix, b::AbstractMatrix) = kron(b, a)

permute(m::AbstractArray, perm) = PermutedDimsArray(m, perm)

function n_qubits_mat_from_external(m::AbstractMatrix)
	(size(m, 1) == size(m, 2)) || error("input matrix is not a square matrix.")
	L = size(m, 1)
	n = round(Int, log2(L))
	(2^n == L) || error("input matrix size error.")
	v = collect(n:-1:1)
	v = vcat(v, v .+ n)
	return reshape(permute(reshape(m, Tuple(2 for i in 1:2*n)), v), size(m))
end

function _kron_ops(op)
	isempty(op) && error("ops is empty.")
	nb = length(op)
	m = op[1]
	for i = 2:length(op)
	    m = kron(m, op[i])
	end
	return m
end

# function contract(a::AbstractArray{Ta, Na}, b::AbstractArray{Tb, Nb}, axs::Tuple{NTuple{N, Int}, NTuple{N, Int}}) where {Ta, Na, Tb, Nb, N}
#     ia, ib = axs
#     seqindex_a = move_selected_index_backward(collect(1:Na), ia)
#     seqindex_b = move_selected_index_forward(collect(1:Nb), ib)
#     ap = permute(a, seqindex_a)
#     bp = permute(b, seqindex_b)
#     return reshape(tie(ap, (Na-N, N)) * tie(bp, (N, Nb-N)), size(ap)[1:(Na-N)]..., size(bp)[(N+1):Nb]...)
# end

# function tie(a::AbstractArray{T, N}, axs::NTuple{N1, Int}) where {T, N, N1}
#     (sum(axs) != N) && error("total number of axes should equal to tensor rank.")
#     return reshape(a, _group_extent(size(a), axs))
# end

# function _group_extent(extent::NTuple{N, Int}, idx::NTuple{N1, Int}) where {N, N1}
#     ext = Vector{Int}(undef, N1)
#     l = 0
#     for i=1:N1
#         ext[i] = prod(extent[(l+1):(l+idx[i])])
#         l += idx[i]
#     end
#     return NTuple{N1, Int}(ext)
# end




# """	
# 	move_selected_index_forward(a, I)
# 	move the indexes specified by I to the front of a
# 	# Arguments
# 	@ a::NTuple{N, Int}: the input tensor.
# 	@ I: tuple or vector of integer.
# """
# function move_selected_index_forward(a::Vector{T}, I) where {T}
#     na = length(a)
#     nI = length(I)
#     b = Vector{T}(undef, na)
#     k1 = 0
#     k2 = nI
#     for i=1:na
#         s = 0
#         while s != nI
#         	if i == I[s+1]
#         		b[s+1] = a[k1+1]
#         	    k1 += 1
#         	    break
#         	end
#         	s += 1
#         end
#         if s == nI
#         	b[k2+1]=a[k1+1]
#         	k1 += 1
#             k2 += 1
#         end
#     end
#     return b
# end

# function move_selected_index_forward(a::NTuple{N, T}, I) where {N, T}
#     return NTuple{N, T}(move_selected_index_forward([a...], I))
# end

# """	
# 	move_selected_index_backward(a, I)
# 	move the indexes specified by I to the back of a
# 	# Arguments
# 	@ a::NTuple{N, Int}: the input tensor.
# 	@ I: tuple or vector of integer.
# """
# function move_selected_index_backward(a::Vector{T}, I) where {T}
# 	na = length(a)
# 	nI = length(I)
# 	nr = na - nI
# 	b = Vector{T}(undef, na)
# 	k1 = 0
# 	k2 = 0
# 	for i = 1:na
# 	    s = 0
# 	    while s != nI
# 	    	if i == I[s+1]
# 	    		b[nr+s+1] = a[k1+1]
# 	    		k1 += 1
# 	    		break
# 	    	end
# 	    	s += 1
# 	    end
# 	    if s == nI
# 	        b[k2+1] = a[k1+1]
# 	        k2 += 1
# 	        k1 += 1
# 	    end
# 	end
# 	return b
# end

# function move_selected_index_backward(a::NTuple{N, T}, I) where {N, T}
# 	return NTuple{N, T}(move_selected_index_backward([a...], I))
# end


