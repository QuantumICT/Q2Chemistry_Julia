


function differentiate(x::AdjointQuantumGate)
	r = differentiate(x.parent)
	isnothing(r) && return r
	# otherwise r should be a list of gates
	return adjoint.(r)
end

# one body parametric gates
function differentiate(x::RxGate)
	if is_parameters(x)[1]
		mm = 0.5 .* Rx(parameters(x)[1] + pi)
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::RyGate)
	if is_parameters(x)[1]
		mm = 0.5 .* Ry(parameters(x)[1] + pi)
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::RzGate)
	if is_parameters(x)[1]
		mm = 0.5 .* Rz(parameters(x)[1] + pi)
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::PHASEGate)
	if is_parameters(x)[1]
		mm = phase_mat_diff(parameters(x)[1])
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
phase_mat_diff(phi::Real) = [0 0; 0. im*exp(im*phi)]

# two body parametric gates
function differentiate(x::CRxGate)
	if is_parameters(x)[1]
		mm = rkron(Gates.DOWN, 0.5 .* Rx(parameters(x)[1] + pi))
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::CRyGate)
	if is_parameters(x)[1]
		mm = rkron(Gates.DOWN, 0.5 .* Ry(parameters(x)[1] + pi))
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::CRzGate)
	if is_parameters(x)[1]
		mm = rkron(Gates.DOWN, 0.5 .* Rz(parameters(x)[1] + pi))
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end
function differentiate(x::CPHASEGate)
	if is_parameters(x)[1]
		mm = rkron(Gates.DOWN, phase_mat_diff(parameters(x)[1]))
		return [QuantumGate(positions(x), mm' * mat(x) )]
	end
	return nothing
end


function differentiate(x::FSIMGate)
	@assert length(parameters(x)) == 5
	paras = parameters(x)
	isparas = is_parameters(x)
	any(isparas) || return nothing
	rs = []
	x_mat = mat(x)
	diffs = [fsim_mat_diff1(paras...), fsim_mat_diff2(paras...), fsim_mat_diff3(paras...), fsim_mat_diff4(paras...), fsim_mat_diff5(paras...)]
	for i in 1:5
		if isparas[i]
			push!(rs, QuantumGate(positions(x), diffs[i]' * x_mat ) )
		end
	end
	return rs
end


function fsim_mat_diff1(theta::Real, phi::Real, deltap::Real, deltam::Real, deltamoff::Real)
	m = [[0 0 0 0]; [0 -exp(im*(deltap+deltam))*sin(theta) -im*exp(im*(deltap-deltamoff))*cos(theta) 0];
    [0 -im*exp(im*(deltap+deltamoff))*cos(theta) -exp(im*(deltap-deltam))*sin(theta) 0]; [0 0 0 0]]
    return n_qubits_mat_from_external(m)		
end

function fsim_mat_diff2(theta::Real, phi::Real, deltap::Real, deltam::Real, deltamoff::Real)
	m = [[0 0 0 0]; [0 0 0 0]; [0 0 0 0]; [0 0 0 -im * exp(im*(2*deltap-phi)) ]]
    return n_qubits_mat_from_external(m)		
end

function fsim_mat_diff3(theta::Real, phi::Real, deltap::Real, deltam::Real, deltamoff::Real)
	m = [[0 0 0 0]; [0 im * exp(im*(deltap+deltam))*cos(theta) exp(im*(deltap-deltamoff))*sin(theta) 0];
    [0 exp(im*(deltap+deltamoff))*sin(theta) im * exp(im*(deltap-deltam))*cos(theta) 0]; [0 0 0 2*im*exp(im*(2*deltap-phi))]]
    return n_qubits_mat_from_external(m)	
end

function fsim_mat_diff4(theta::Real, phi::Real, deltap::Real, deltam::Real, deltamoff::Real)
	m = [[0 0 0 0]; [0 im * exp(im*(deltap+deltam))*cos(theta) 0 0];
    [0 0 -im*exp(im*(deltap-deltam))*cos(theta) 0]; [0 0 0 0]]
    return n_qubits_mat_from_external(m)	
end

function fsim_mat_diff5(theta::Real, phi::Real, deltap::Real, deltam::Real, deltamoff::Real)
	m = [[0 0 0 0]; [0 0 -exp(im*(deltap-deltamoff))*sin(theta) 0];
    [0 exp(im*(deltap+deltamoff))*sin(theta) 0 0]; [0 0 0 0]]
    return n_qubits_mat_from_external(m)	
end





