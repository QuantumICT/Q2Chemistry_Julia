# gradient of expectation values


function apply!(m::QubitsTerm, state::MPS; trunc::TruncationScheme=DefaultMPSTruncation)
	errs = [QuantumSpins._apply_impl((k,), Array(v), state, trunc) for (k, v) in zip(QuantumCircuits.positions(m), oplist(m))]
	state[1] *= QuantumCircuits.coeff(m)
	return errs
end

@adjoint expectation(m::QubitsTerm, state::MPS; kwargs...) = _qterm_expec_util(m, state; kwargs...)
@adjoint expectation(m::QubitsOperator, state::MPS; kwargs...) = _qop_expec_util(m, state; kwargs...)

function _qterm_expec_util(m::QubitsTerm, state::MPS; trunc::TruncationScheme=DefaultMPSTruncation)
	return expectation(m, state), z -> begin
		if ishermitian(m)
			r = copy(state)
			apply!((2 * real(z)) * m, r; trunc=trunc)
			canonicalize!(r, normalize=false, alg=SVDFact(trunc=trunc))
		else
			state_a = copy(state)
			state_b = copy(state)
			apply!(conj(z) * m, state_a; trunc=trunc)
			apply!(z * m', state_b; trunc=trunc)
			# r = state_a + state_b
			alg = StableArith(D=trunc.D, ϵ=trunc.ϵ)
			r, err = mpsadd([state_a, state_b], alg)
		end
		return (nothing, r)
	end 
end

function _qop_expec_util(m::QubitsOperator, state::MPS; trunc::TruncationScheme=DefaultMPSTruncation)
	return expectation(m, state), z -> begin
		alg = StableArith(D=trunc.D, ϵ=trunc.ϵ)
		_is_herm = true #ishermitian(m)
		mpo = _MPO(length(state), m)
		if _is_herm
			mpsout, err = mpompsmult(mpo, state, alg)
			mpsout[1] .*= (2 * real(z))
		else
			a, err = mpompsmult(mpo, state, alg)
			a[1] .*= conj(z)
			b, err = mpompsmult(mpo', state, alg)
			b[1] .*= z
			mpsout, err = mpsadd([a, b], alg)
		end
		return (nothing, mpsout)
	end 
end

function _prodmpo(physpaces::Vector{Int}, m::QubitsTerm) 
	return prodmpo(eltype(m), physpaces, QuantumCircuits.positions(m), oplist(m)) * QuantumCircuits.coeff(m)
end

function _qterms(x::QubitsOperator) 
	r = []
	for (k, v) in x.data
		for (m, c) in v
			if c != zero(c)
				a = QubitsTerm(k, m, c)
				push!(r, a)
			end
		end
	end
	return r
end

function _MPO(L::Int, h::QubitsOperator; alg::MPOCompression=Deparallelise())
	physpaces = [2 for i in 1:L]
	local mpo
	compress_threshold = 20
	for m in _qterms(h)
		if @isdefined mpo
			mpo += _prodmpo(physpaces, m)
		else
			mpo = _prodmpo(physpaces, m)
		end
		if bond_dimension(mpo) >= compress_threshold
			mpo = compress!(mpo, alg=alg)
			compress_threshold += 5
		end
	end
	mpo = compress!(mpo, alg=alg)
	return mpo
end
