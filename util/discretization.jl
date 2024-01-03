using LinearAlgebra


function vanLoan(F::Matrix, G::Matrix; dt::Number=1.0)
    
    n = size(F, 1)
    A = zeros(2*n, 2*n)
    
    A[1:n, 1:n] = -F * dt
    A[1:n, n+1:end] = G * G' * dt
    A[n+1:end, n+1:end] = F' * dt
    
    B = exp(A)
    sigma = B[n+1:end, n+1:end]'
    Q = sigma * B[1:n, n+1:end]
    Q = .5 * (Q + Q')
    return sigma, Q
end

function analyticQ(Mi,λ,σ²; Δt=1)
    """
    Analyic derivation of ∫_0^Δt F(t)QcF(t)'dt
    where F(t) = I + At and Qc = [0,1]*qc*[0 1]
    """
    dims = size(Mi,1)

    # Spectral noise of Matern-1/2 kernels
    qc = 2*λ*σ²

    # Blocks of process noise covariance matrix
    B1 =  Δt^3/3.0 *qc*Mi
    B2 = (Δt^2/2.0 - λ*Δt^3/3.0) *qc*Mi
    B3 = qc*Δt*(1 - λ*Δt + λ^2*Δt^2 ./3.0)*diagm(ones(dims))
    
    return [B1 B2; B2' B3]
end