###
### DQGMRES in julia
###


function incompleteArnoldi(v0::Array{Float64,1}, A::Array{Float64,2}, k::Int, m::Int)

  V = Array{Float64,2}(length(v0),m)
  H = Array{Float64,2}(length(v0),m)

  V[:,1] = v0
  for j = 1:m

    w = A*V[:,j]
    for i = max(1,j-k+1):j
    
      h[i,j] = w'*V[:,i]
      w -= H[i,j]V[:,i]

    end
    h[j+1,j] = norm(w)
    v[:,j+1] = w/H[j+1,j]

  end

  return [H[:,], V[:,m+1]]

end


function dqgmres(A::Array{Float64,2},
	             b::Array{Float64,1},
	             k::Int;
	             tol::Float64=1e-6,
	             maxit::Int=min(length(b),20),
	             x0::Array{Float64,1}=zeros(length(b)))

  γ = Array{Float64,1}(maxit+1)
  V = Array{Float64,2}(length(b), maxit+1)
  H = zeros(maxit+1, maxit+1)
  c = Array{Float64,1}(maxit+1)
  s = Array{Float64,1}(maxit+1)
  P = Array{Float64,2}(length(b), maxit+1)

  x      = x0
  r      = b - A*x
  γ[1]   = norm(r)
  V[:,1] = r/γ[1]

  for m = 1:maxit

    l = max(1,m-k+1)
    w = A*V[:,m]
    H[1:m,m] = w'*V[:,1:m]
    w -= V[:,1:m]*H[1:m,m]
    # for i = l:m

    #   H[i,m] = w'*V[:,i]
    #   w      -= H[i,m]*V[:,i]

    # end
    H[m+1,m] = norm(w)
    V[:,m+1] = w/H[m+1,m]

    for i = l:m-1

      H[i,m] = c[i]*H[i,m] + s[i]*H[i+1,m]

    end

    c[m] = H[m,m]/sqrt(H[m,m]^2 + H[m+1,m]^2)
    s[m] = H[m+1,m]/sqrt(H[m,m]^2 + H[m+1,m]^2)

    γ[m+1] = -s[m]*γ[m]
    γ[m]   = c[m]*γ[m]
    H[m,m] = c[m]*H[m,m] + s[m]*H[m+1,m]
    H[m+1,m] = 0

    # P[:,m] = (V[:,m] - P[:,l:m-1]*H[l:m-1,m])/H[m,m]
    # P[:,m] = (V[:,m] - Float64.(sum([H[i,m]*P[:,i] for i in l:m-1])))/H[m,m]    

    # x = x + γ[m]*P[:,m]

    # if abs(γ[m+1]) < tol

      # return V, H, c, s, P, γ

    # end

  end

  # return x, γ[m+1], maxit
  V, H, c, s, P, γ, x

end