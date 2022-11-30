
function expectation(psiA::MPS, m::QubitsTerm, psiB::MPS, envs::OverlapCache=environments(psiA, psiB))
	cstorage = envs.cstorage
	(length(psiA) == length(psiB) == length(cstorage)-1) || throw(DimensionMismatch())
	pos = QuantumCircuits.positions(m)
	isempty(pos) && return 0.
	L = length(psiA)
	ops = oplist(m)
	pos_end = pos[end]
	(pos_end <= L) || throw(BoundsError())
	@tensor hold[-1; -3] := conj(psiA[pos_end][-1, 1, 2]) * cstorage[pos_end+1][2, 3] * psiB[pos_end][-3, 5, 3] * ops[end][1, 5]
	for j in pos_end-1:-1:1
		pj = findfirst(x->x==j, pos)
		if isnothing(pj)
			hold = updateright(hold, psiA[j], pj, psiB[j])
		else
			hold = updateright(hold, psiA[j], ops[pj], psiB[j])
		end
	end
	return QuantumSpins.scalar(hold) * QuantumCircuits.coeff(m)
end

expectation(m::QubitsTerm, psi::MPS, envs::OverlapCache=environments(psi, psi); kwargs...) = expectation(psi, m, psi, envs)

function expectation(psiA::MPS, h::QubitsOperator, psiB::MPS)
	(QuantumCircuits.get_largest_pos(h) <= length(psiA)) || throw(DimensionMismatch())
	envs = environments(psiA, psiB)
	r = zero(promote_type(scalar_type(psiA), scalar_type(psiB)))
	for (k, v) in h.data
		for (x, c) in v
			m = QubitsTerm(k, x, c)
			r += expectation(psiA, m, psiB, envs)
		end
	end
	return r
end

expectation(h::QubitsOperator, psi::MPS; kwargs...) = expectation(psi, h, psi)


# density operator
"""
	expectation(m::QubitsTerm, psi::DensityOperatorMPS) 
	⟨I|h|ρ⟩
"""
function expectation(m::QubitsTerm, psi::DensityOperatorMPS, envs=environments(psi))
	sm = QubitsTerm(QuantumCircuits.positions(m), [QuantumSpins.rkron(a, one(a)) for a in oplist(m)], QuantumCircuits.coeff(m))
	return expectation(psi.I, sm, psi.data, envs)
end

function expectation(h::QubitsOperator, psi::DensityOperatorMPS)
	(length(h) <= length(psi)) || throw(DimensionMismatch())
	envs = environments(psi, psi)
	r = zero(scalar_type(psi))
	for (k, v) in h.data
		for (x, c) in v
			m = QubitsTerm(k, x, c)
			r += expectation(m, psi, envs)
		end
	end
	return r
end


