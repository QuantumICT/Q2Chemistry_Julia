

# density matrix initializer
function densitymatrix(state::MPS; trunc::TruncationScheme=DefaultMPSTruncation) 
	@assert all(physical_dimensions(state) .== 2)
	return DensityOperator(state; trunc=trunc)
end

densitymatrix_mps(::Type{T}, states::Vector{Int}; kwargs...) where {T <: Number} = densitymatrix(statevector_mps(T, states); kwargs...)
densitymatrix_mps(::Type{T}, n) where {T <: Number} = densitymatrix_mps(T, zeros(Int, n))
densitymatrix_mps(n::Int) = densitymatrix_mps(ComplexF64, n)

