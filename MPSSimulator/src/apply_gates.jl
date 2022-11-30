
function apply!(s::Gate, mps::MPS; trunc::TruncationScheme=DefaultMPSTruncation) 
	(length(QuantumCircuits.positions(s)) <= 4) || throw(ArgumentError("only 4-body (or less) gates are currently allowed."))
	QuantumSpins.svectors_uninitialized(mps) && canonicalize!(mps, normalize=false, alg=SVDFact(trunc=trunc))
	return [QuantumSpins._apply_impl(QuantumCircuits.ordered_positions(s), Array(QuantumCircuits.ordered_op(s)), mps, trunc)]
end 

function apply!(s::Gate, rho::DensityOperatorMPS; trunc::TruncationScheme=DefaultMPSTruncation)
	(length(QuantumCircuits.positions(s)) <= 4) || throw(ArgumentError("only 4-body (or less) gates are currently allowed."))
	QuantumSpins.svectors_uninitialized(rho.data) && canonicalize!(rho, normalize=false, alg=SVDFact(trunc=trunc))
	mop = QuantumCircuits.ordered_op(s)
	mop = QuantumSpins.rkron(mop, conj(mop))
	return [QuantumSpins._apply_impl(QuantumCircuits.ordered_positions(s), mop, rho.data, trunc)]
end

function apply!(circuit::QCircuit, mps::AbstractMPS; kwargs...)
	errs = Float64[]
	for g in circuit
		err = apply!(g, mps; kwargs...)
		append!(errs, err)
	end
	return errs
end


function apply(circuit::QCircuit, mps::AbstractMPS; kwargs...) 
	mpsc = copy(mps)
	apply!(circuit, mpsc; kwargs...)
	return mpsc
end