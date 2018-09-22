#=
 julia adaptation of MATLAB's GMRES implementation
 -------------------------------------------------

 it works now! but it's slow. consider:
 	- overall worksace ~ ~ ~ struct probably doesn't need inputs
=#
using LinearAlgebra
using SparseArrays

# put it all in a big ol' struct
mutable struct GMRESIterable{T, matT}
	# algorithmic workspace
	J::Matrix{T}	# givens rotation coefficients
	U::Matrix{T}	# householder reflector vectors
	R::Matrix{T}	# krylov coefficients (?) [[this is the trick that gives us cheap updates]]
	w::Vector{T}	# residiual norm tracker
    k::Int			# inner loop index

	# computational workspace
	A::matT
    Ml::matT		# left preconditioner
    Mr::matT		# right preconditioner
    b::Vector{T}		# right-hand-side

    restart::Int
    maxiter::Int
    tol::T

	x::Vector{T}		# solution
    u::Vector{T}		# working vector
    r::Vector{T}		# residual vector

    evalxm::Bool		# track whether solution has been updated
    moresteps::Bool		# moar steps

    resvec::Vector{T}	# residual norm history
end

done(g::GMRESIterable, iteration::Int) =  iteration ≥ g.maxiter || g.tol ≥ norm(g.r)

# loops appear here
function Base.iterate(g::GMRESIterable, state=0)
	# @printf "pls"
	iteration = state
	if done(g, iteration)
		# @printf "done!"
		return nothing
	end

	# @show g.k
	# @show g.u[1:11]

	v = -2g.u[g.k] * g.u
	v[g.k] += 1

	# @show g.k
	# @show v[1:11]

	for j = g.k-1:-1:1
		uu = deepcopy(g.U[:,j])
		v -= 2uu'*v * uu
	end

	# explicitly normalize for accuracy before applying A
	normalize!(v)
	v = g.A * v

	# ***preconditioner stuff goes here

	# form Pj*Pj-1*...*P1*A*v
	for j = 1:g.k
		uu = deepcopy(g.U[:,j])
		v -= 2uu'*v * uu
	end

	# determine Pj+1
	if g.k != g.restart+1
		g.u = deepcopy(v)
		g.u[1:g.k] = zero(g.u[1:g.k])
		# @show g.k
		α = norm(g.u)
		if α != 0.
			α *= sign(v[g.k+1])
			g.u[g.k+1] += α
			normalize!(g.u)
			if g.k != g.restart
				g.U[:,g.k+1] = deepcopy(g.u)
			end

			# apply Pj+1 to v
			v[g.k+2:end] = zero(v[g.k+2:end])
			v[g.k+1] = -α
		end
	end

	# @show g.k
	# @show g.u[1:11]

	# apply givens rotations to v
	for j = 1:g.k-1
		vv = deepcopy(v[j])
		v[j] = conj(g.J[1,j])*v[j] + conj(g.J[2,j])*v[j+1]
		v[j+1] = -g.J[2,j]*vv + g.J[1,j]*v[j+1]
	end

	# @show g.k
	# @show v[1:11]

	# compute givens rotation Jm
	if g.k != g.restart+1
		# @show g.k
		ρ = norm(v[g.k:g.k+1])
		g.J[:,g.k] = v[g.k:g.k+1] / ρ
		g.w[g.k+1] = -g.J[2,g.k] * g.w[g.k]
		g.w[g.k] *= conj(g.J[1,g.k])
		v[g.k] = deepcopy(ρ)
		v[g.k+1] = 0
	end

	g.R[:,g.k] = deepcopy(v[1:g.restart])

	normr = abs(g.w[g.k+1])
	g.resvec[iteration+2] = deepcopy(normr)

	g.k += 1

	###***###***###***###*** end inner loop

	# compute x at restart/end
	if g.k == g.restart+1
		g.k -= 1

		if !g.evalxm
			y = g.R[1:g.k,1:g.k] \ g.w[1:g.k]
			additive = g.U[:,g.k] * -2y[g.k]*conj(g.U[g.k,g.k])
			additive[g.k] += y[g.k]

			for j = g.k-1:-1:1
				additive[j] += y[j]
				additive -= g.U[:,j] * 2g.U[:,j]'*additive
			end

			g.x += additive
			g.evalxm = true
		else
			addvc = [ -(g.R[1:g.k-1,1:g.k-1] \ g.R[1:g.k-1,g.k]) * g.w[g.k] / g.R[g.k,g.k];
					   g.w[g.k] / g.R[g.k,g.k] ]
			additive = g.U[:,g.k] * -2addvc[g.k]*conj(g.U[g.k,g.k])
			additive[g.k] += addvc[g.k]
			for j = g.k-1:-1:1
				additive[j] += addvc[j]
				additive -= g.U[:,j] * 2g.U[:,j]'*additive
			end
			g.x += additive
		end

		g.r = g.b - g.A*g.x
		g.resvec[iteration+1] = norm(g.r)

		# restart when not done
		if !done(g, iteration)
			g.J = zero(g.J)
			g.U = zero(g.U)
			g.R = zero(g.R)
			g.w = zero(g.w)
			g.k = 1

			g.u = deepcopy(g.r)
			normr = norm(g.r)
			β = sign(g.r[1]) * normr
			g.u[1] += β
			normalize!(g.u)

			g.U[:,1] = deepcopy(g.u)
			g.w[1] = -β

			g.evalxm = false
			# rv = deepcopy(g.resvec)
			# g = GMRESIterable!(g.x, g.A, g.b;
			# 				   tol=g.tol,
			# 				   restart=g.restart,
			# 				   maxiter=g.maxiter)
			# g.resvec = rv
		end
		# @show g.k
	elseif g.resvec[iteration+1] <= g.tol #* norm(g.b)
		if !g.evalxm
			y = g.R[1:g.k,1:g.k] \ g.w[1:g.k]
			additive = g.U[:,g.k] * -2y[g.k]*conj(g.U[g.k,g.k])
			additive[g.k] += y[g.k]
			for j = g.k-1:-1:1
				additive[j] += y[j]
				additive -= g.U[:,j] * 2g.U[:,j]'*additive
			end
			g.x += additive
			g.evalxm = true
		else
			addvc = [ -(g.R[1:g.k-1,1:g.k-1] \ g.R[1:g.k-1,g.k]) * g.w[g.k] / g.R[g.k,g.k];
					   g.w[g.k] / g.R[g.k,g.k] ]
			additive = g.U[:,g.k] * -2addvc[g.k]*conj(g.U[g.k,g.k])
			additive[g.k] += addvc[g.k]
			for j=g.k-1:-1:1
				additive[j] += addvc[j]
				additive -= g.U[:,j] * 2g.U[:,j]'*additive
			end
			g.x += additive
		end
		g.r = g.b - g.A*g.x
	end
	# @show g.r[1:10]

	g.x, iteration + 1
end

# iterable object constructor
function GMRESIterable!(x, A, b;
						Ml=sparse(I,size(A)), Mr=sparse(I,size(A)),
						tol=1e-6,
						restart::Int=min(20,size(A,1)),
						maxiter::Int=size(A,1)) where {matT}
	T = eltype(A)
	J = zeros(T, 2, restart)
	U = zeros(T, size(A,1), restart)
	R = zeros(T, restart, restart)
	w = zeros(T, restart + 1)
	resvec = Inf * ones(T, maxiter + 1)

	r = b - A*x
	u = deepcopy(r)
	normr = norm(r)
	β = sign(r[1])*normr
	resvec[1] = deepcopy(β)
	u[1] += β
	normalize!(u)

	U[:,1] = u
	w[1] = -β

	GMRESIterable{T,typeof(A)}(J, U, R, w, 1,

						  A, Ml, Mr, b,

						  restart, maxiter, tol,

						  x, u, r,

						  false, false,

						  resvec)
end

# function wrappers
