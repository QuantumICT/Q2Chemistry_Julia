


@adjoint apply(circuit::QCircuit, x::AbstractMPS; kwargs...) = begin
    y = apply(circuit, x; kwargs...)
    return y, Δ -> begin
        Δ, grads, y = back_propagate(copy(Δ), circuit, copy(y); kwargs...)
        return grads, Δ
    end
end

@adjoint qubit_encoding_mps(::Type{T}, mpsstr::Vector{<:Real}) where {T<:Number} = begin
    y = qubit_encoding_mps(T, mpsstr)
    return y, Δ -> begin
        circuit = QCircuit([RyGate(i, theta*pi, isparas=true) for (i, theta) in enumerate(mpsstr)])
        Δ, grads, y = back_propagate(Δ, circuit, copy(y))
        return nothing, grads .* pi
    end
end


function _gate_expec(x::AbstractMPS, m::Gate, y::AbstractMPS; kwargs...)
    yc = copy(y)
    apply!(m, yc; kwargs...)
    return dot(x, yc)
end

function back_propagate(Δ::AbstractMPS, m::Gate, y::AbstractMPS; kwargs...)
    Δerr = apply!(m', Δ; kwargs...)
    yerr = apply!(m', y; kwargs...)
    ∇θs = nothing
    if nparameters(m) > 0
        ∇θs = [real(_gate_expec(y, item, Δ; kwargs...)) for item in differentiate(m)]
    end
    return Δ, ∇θs, y
end


function back_propagate(Δ::AbstractMPS, circuit::QCircuit, y::AbstractMPS; kwargs...)
    RT = real(scalar_type(y))
    grads = Vector{RT}[]
    for item in reverse(circuit)
        Δ, ∇θs, y = back_propagate(Δ, item, y; kwargs...)
        !isnothing(∇θs) && push!(grads, ∇θs)
    end

    ∇θs_all = RT[]
    for item in Iterators.reverse(grads)
        append!(∇θs_all, item)
    end

    return Δ, ∇θs_all, y
end


# function back_propagate(Δ::AbstractMatrix, m::Gate, y::DensityMatrix)
#     Δ = StateVector(Δ, nqubits(y))
#     Δ = apply!(m', Δ)
#     y = apply!(m', y)
#     ∇θs = nothing
#     if nparameters(m) > 0
#         ∇θs = [real(expectation(y, item, Δ)) for item in differentiate(m)]
#     end
#     return storage(Δ), ∇θs, y
# end