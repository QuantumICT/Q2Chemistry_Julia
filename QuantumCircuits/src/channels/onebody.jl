

# predefined quantum channels

# see NielsenChuang
AmplitudeDamping(pos::Int; γ::Real) = QuantumMap(1, [[1 0; 0 sqrt(1-γ)], [0 sqrt(γ); 0 0]])
PhaseDamping(pos::Int; γ::Real) = QuantumMap(1, [[1 0; 0 sqrt(1-γ)], [0 0; 0 sqrt(γ)]])
Depolarizing(pos::Int; p::Real) = QuantumMap(1, [sqrt(1-3*p/4) .* I₂, (sqrt(p)/2) .* X, (sqrt(p)/2) .* Y, (sqrt(p)/2) .* Z ])