function [ omega, obj ] = SPCA_PSD_modified( S , omega, lambda, eta)

[ d , ~] = size(S);
Id = eye(d);

delta = inf;
obj = zeros(1,50);
k = 1;
obj(k) = abs(trace(S-omega*(S+S')+omega*S*omega')) + lambda*sum((sqrt(sum(abs(omega).^2)))) + eta * trace(omega);
Niter = 50;

while delta > 10^-5
    diag_W = 0.5 * (sqrt(sum(abs(omega).^2)+eps)).^(-1);
    W = diag(diag_W);
    
    omega = conj(( (S+S')/2 - eta / 2 * Id)/( S + lambda * W + 0.001*eye(d) ));
    omega = Keep_PSD(omega);
    
    obj(k+1) = abs(trace(S-omega*(S+S')+omega*S*omega')) + lambda*sum((sqrt(sum(abs(omega).^2)))) + eta * trace(omega);
    delta = abs(obj(k+1)-obj(k)); 
    k = k + 1;
    
    if k > Niter
        break
    end

end

end




