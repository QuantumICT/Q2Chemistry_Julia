
abstract type ParametricGate{N} <: Gate{N} end


parameters(x::Gate) = nothing
is_parameters(x::Gate) = nothing
active_parameters(x::Gate) = nothing
nparameters(x::Gate) = 0

parameters(x::ParametricGate) = x.paras
parameters(x::AdjointQuantumGate) = parameters(x.parent)
is_parameters(x::ParametricGate) = x.isparas
is_parameters(x::AdjointQuantumGate) = is_parameters(x.parent)
nparameters(x::ParametricGate) = sum(is_parameters(x))
active_parameters(x::ParametricGate) = x.paras[x.isparas]
active_parameters(x::AdjointQuantumGate) = active_parameters(x.parent)
nparameters(x::AdjointQuantumGate) = nparameters(x.parent)
activate_parameter!(x::Gate, j::Int) = x
activate_parameters!(x::Gate) = x
deactivate_parameter!(x::Gate, j::Int) = x
deactivate_parameters!(x::Gate) = x

function reset_parameters_util!(x::Gate, p::Vector{<:Real}, pos::Int)
	n = nparameters(x)
	(n == 0) && return pos
	parameters(x)[is_parameters(x)] .= view(p, pos+1:(pos+n))
	return pos + n
end

function activate_parameter!(x::ParametricGate, j::Int)
	x.isparas[j] = true
	return x
end
function activate_parameters!(x::ParametricGate)
	x.isparas .= true
	return x
end
function deactivate_parameter!(x::ParametricGate, j::Int)
	x.isparas[j] = false
	return x
end
function deactivate_parameters!(x::ParametricGate)
	x.isparas .= false
	return x
end

function activate_parameter!(x::AdjointQuantumGate, j::Int)
	activate_parameter!(x.parent, j)
	return x
end
function activate_parameters!(x::AdjointQuantumGate)
	activate_parameters!(x.parent)
	return x
end
function deactivate_parameter!(x::AdjointQuantumGate, j::Int)
	deactivate_parameter!(x.parent, j)
	return x
end
function deactivate_parameters!(x::AdjointQuantumGate)
	deactivate_parameters!(x.parent)
	return x
end

differentiate(x::Gate) = nothing

"""
	differentiate(x::ParametricGate)
	return a list of gates, with the same number as the active parameters
"""
differentiate(x::ParametricGate) = error("differentiate not implemented for parametric gate $(typeof(x)).") 

