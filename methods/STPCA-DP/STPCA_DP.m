function [output] = STPCA_DP(X,para,para_special)
% Sparse Tensor PCA based on direction-unfolding product
% X: a data tensor of size fea1*fea2*..*Sample
% S: a cell that contains directions along which the data matrix is
%    reconstructed, each direction is represented as a row vector. 
% lamda: a vector that contains each regularization parameter of l2,1-norm
%        for each direction
% eta: a vector that contains each regularization parameter of nuclear norm
%      for each direction
lambda = para.lambda;
eta = para.eta;
S = para_special.S;
%%%% examining input parameters
K1 = size(S,1);
K2 = size(S,2);
if (K1~=1 && K2~=1) || (iscell(S) == 0)
    error('S should be a 1*K cell. ');
end
if length(unique((cell2mat(S))))~=length((cell2mat(S)))
   error('Repeated directions occurred') 
end
K = length(S);
lambda = repmat(lambda,[1,K]);
eta = repmat(eta,[1,K]);



%%%% Initialization
size_X = size(X);
nSamp = size_X(end);
X_mean = repmat(mean(X,length(size_X)),[ones(1,length(size_X)-1),nSamp]);
X_centralized = X-X_mean; % Centralization
X_ten = tensor(X_centralized);
X_RecS = X_ten;
num_modes = ndims(X)-1;
Omega = cell(K,1);
size_omega_ten = cell(K,1);
l21_norm = zeros(K,1);
nuclear_norm = zeros(K,1);

for k = 1:K
    omega = diag(rand(prod(size_X(S{k})),1));
    Omega{k} = tensor(omega,[prod(size_X(S{k})),prod(size_X(S{k}))]); % Intialization of each mode's omega
    size_omega_ten{k} = size(Omega{k});
    X_RecS_mode = tenmat(X_RecS,S{k});
    X_RecS_mode(:,:) = tenmat(Omega{k},1:(length(size_omega_ten{k})/2))*X_RecS_mode; % S-direction product
    X_RecS = tensor(X_RecS_mode);
    l21_norm(k) = sum((sqrt(sum(abs(omega).^2))));
    nuclear_norm(k) = trace(omega);
end
flag = 1;
fros = abs(double(X_ten-X_RecS)).^2;
OBJ(flag) = sum(fros(:)) + lambda*l21_norm + eta*nuclear_norm;

%%%% Optimization

delta = inf;
while delta >10^-5  
    X_RecS = X_ten;
    for k = 1:K
        X_RecS_other = X_ten;
        for kk = 1:K
           if kk ~=k
               X_RecS_mode_other = tenmat(X_RecS_other,S{kk});
               X_RecS_mode_other(:,:) = tenmat(Omega{kk},1:(length(size_omega_ten{kk})/2))*X_RecS_mode_other;
               X_RecS_other = tensor(X_RecS_mode_other);
           end
        end
        omega = double(tenmat(Omega{k},1:(length(size_omega_ten{k})/2)));
        X_RecS_other_mode = double(tenmat(X_RecS_other, S{k}));
        X_ten_mode = double(tenmat(X_ten, S{k}));
        Cov_mode = X_ten_mode*X_RecS_other_mode'; % sum the covariance matrices of all samples
        [ omega, ~ ] = SPCA_PSD_modified( Cov_mode , omega, lambda(k), eta(k));
        Omega{k} = tensor(omega,size_omega_ten{k});
        X_RecS_mode = tenmat(X_RecS,S{k});
        X_RecS_mode(:,:) = tenmat(Omega{k},1:(length(size_omega_ten{k})/2))*X_RecS_mode;
        X_RecS = tensor(X_RecS_mode);
        l21_norm(k) = sum((sqrt(sum(abs(omega).^2))));
        nuclear_norm(k) = trace(omega);
    end
    
    flag = flag + 1;
    fros = abs(double(X_ten-X_RecS)).^2;
    OBJ(flag) = sum(fros(:)) + lambda*l21_norm + eta*nuclear_norm;
    delta = abs(OBJ(flag)-OBJ(flag-1));
    if flag > 100 || (K==1)
        break
    end
end

%%%% Scoring

OMEGA = double(Omega{end});
score_EachMode = cell(K,1);
id_EachMode = cell(K,1);
for k = 1:K
   score_EachMode{k} = sum(abs(double(Omega{k})).^2);
   [~,id_EachMode{k}] = sort(score_EachMode{k},'descend');
end
if K > 1
    for k = K-1:-1:1
        OMEGA = kron(OMEGA,double(Omega{k}));
    end
end
score = (sum(abs(OMEGA).^2))';
if K == 1
   rest_size = size_X;
   rest_size(S{1}) = []; 
   rest_size(end) = [];
   score = repmat(score,[prod(rest_size),1]);
end

[~,id] = sort(score,'descend');
output.id = id;
output.score = score;
output.id_EachMode = id_EachMode;
output.score_EachMode = score_EachMode;
output.obj = OBJ;
end

