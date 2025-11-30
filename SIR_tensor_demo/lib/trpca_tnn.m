function [L,S,err,iter] = trpca_tnn(X,lambda,mu,max_iter,tol)

% Solve the Tensor Robust Principal Component Analysis based on Tensor Nuclear Norm problem by ADMM
%
% min_{L,S} ||L||_*+lambda*||S||_1, s.t. X=L+S
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations
%
% version 1.0 - 19/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

rho = 1.1;
max_mu = 1e10;
DEBUG = 0;


dim = size(X);
L = zeros(dim);
S = L;
Y = L;

iter = 0;
for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % update L
    [L,tnnL] = prox_tnn(-S+X-Y/mu,1/mu);

    
    % update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);
  
    dY = L+S-X;
    
    chg = norm(L(:)-Lk(:))/norm(Lk(:));
   
        if iter == 10 || mod(iter, 10) == 0
              disp(['at iteration ',num2str(iter), '\ chg ',num2str(chg)]) 
        end

    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
end
err = norm(dY(:));
