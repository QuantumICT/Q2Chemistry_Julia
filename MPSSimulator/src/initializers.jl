

# initializer
statevector_mps(::Type{T}, states::Vector{Int}) where {T<:Number} = prodmps(T, [2 for i in 1:length(states)], states)
statevector_mps(::Type{T}, n::Int) where {T <: Number} = statevector_mps(T, zeros(Int, n))
statevector_mps(n::Int) = statevector_mps(ComplexF64, n)

kernal_mapping(s::Real) = [cos(s*pi/2), sin(s*pi/2)]

qubit_encoding_mps(::Type{T}, v::AbstractVector{<:Real}) where {T<:Number} = prodmps(T, kernal_mapping.(v))
qubit_encoding_mps(v::AbstractVector{<:Real}) = qubit_encoding_mps(ComplexF64, v)
