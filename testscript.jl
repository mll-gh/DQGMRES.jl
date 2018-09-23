###
### numerical testing script
###

using DataFrames
using CSV
using PyPlot

# for reading COO-sparse format csv files
@inline function csvsparse(file)
	A_ = CSV.read(file)
	n = size(A_,1)
	A = spzeros(n,n)
	for z in 1:n
		A[ A_[z,1], A_[z,2] ] = A_[z,3]
	end
	A
end

# test method convergence on sparse matrix A for a fixed tol & k
function conv_plot(A::AbstractMatrix, k::Int, tol::Float64, maxiter=Int)
	# setup for solvers
	b = ones(size(A,1))
	# return residual vectors
	rsv_gm = gmres(A,b; restart=k, tol=tol, maxiter=maxiter)
	# rsv_dq = dqgmres(A,b; k=k, tol=tol, maxiter=maxiter)
	## generate plot
	plot(rsv_gm, label="GMRES")
	# plot(rsv,dq, label="DQGMRES")
	xlabel("iterations")
	ylabel("residual norm")
	title("convergence [k=$k]")
	legend()
end

function stag(rvc)
	n = length(rvc)
	s = 0
	for z in 1:n-1
		s = rvc[z+1] == rvc[z] ? s+1 : 1
		if s > 2
			return z-1 #, rvc[z-1]
		end
	end
	n
end

# overall parameter scan
function big_conv_plot(A::AbstractMatrix, tol::Float64, maxiter::Int)
	# solver/plot setup
	n = size(A,1)
	b = ones(n)
	steps_gm = zeros(n)
	steps_dq = zeros(n)
	for k in 1:n
		rsv_gm = gmres(A,b; restart=k, tol=tol, maxiter=maxiter)
		# rsv_dq = dqgmres(A,b; k=k, tol=tol, maxiter=maxiter)
		steps_gm[k] = stag(rsv_gm)
		# steps_dq[k] = stag(rsv_dq)
	end
	## generate plot
	scatter(1:n, steps_gm, label="GMRES")
	# plot(steps_dq, label="DQGMRES")
	xlabel("k")
	ylabel("steps")
	title("convergence [tol=$tol]")
	legend()
	steps_gm
end
