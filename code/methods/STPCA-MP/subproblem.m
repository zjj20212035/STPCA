function [omega, id, sumomega, obj, time] = subproblem( X, lambda, eta)
% Solving each subproblem in STPCA-MP

% X: data matrix, each column is a sample
% lambda: the regularization parameter for l2,1-norm
% eta: the regularization parameter for trace function
% omega: The reconstruction matrix
% id: The rank of features
% obj: The objective function value

[ d , n ] = size( X );
Id = eye(d);
In = eye(n);
X = X - repmat(mean(X,2),1,n);
S = X*X';
X_H = X';

delta = inf;
omega = diag(rand(1,d));
obj = zeros(1,50);
k = 1;
obj(k) = norm(X-omega*X,'fro')^2 + lambda*sum((sqrt(sum(abs(omega).^2)))) + eta * trace(omega);
Niter = 100;
time = [];

while delta > 10^-5
    tic
    diag_W = 0.5 * (sqrt(sum(abs(omega).^2)+eps)).^(-1);
    W = diag(diag_W);
    iD = diag(1./(lambda*diag_W + 0.01));
    
    if d < n
        omega = conj(( S - eta / 2 * Id)/( S + lambda * W + 0.001*eye(d) ));
    elseif d >= n
        omega = iD * X / (In + X_H *iD*X) * X_H .* repmat(diag(Id + eta/2*iD)',d,1) - eta/2*iD;
    end
    omega = Keep_HPSD(omega);
    
    obj(k+1) = norm(X-omega*X,'fro')^2 + lambda*sum((sqrt(sum(abs(omega).^2)))) + eta * trace(omega);
    delta = abs(obj(k+1)-obj(k)); 
    k = k + 1;
    
    if k > Niter
        break
    end
    time = [time,toc];

end

sqomega = (abs(omega).^2);
sumomega = sum(sqomega);
[~,id] = sort(sumomega,'descend');
end


