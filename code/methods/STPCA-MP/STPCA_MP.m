function [output] = STPCA_MP(X_ten,para,para_special)
% Sparse Tensor PCA based on *M product
% X_ten (d1*d2*nSamp): a data tensor that contains a series of 2D samples
% M: an invertible matrix, or "fft" (denotes Fast Fourier Transform) and
%   "I" (denotes identity matrix)
% lambda, eta: regularization parameters
% orientation: specifying order set (1: D^o={1,3,2} or 2: D^o={2,3,1})
lambda = para.lambda;
eta = para.eta;
tensor_size = para_special.tensor_size;
M = para_special.M;
orientation = para_special.orientation;
if M~="fft" && M~="I"
   M_matrix = para_special.M_matrix;
end

nFea = prod(tensor_size);
size_X = size(X_ten);

if orientation == 1
    %% if direction of 1-mode is chosen
    X_ten = permute(X_ten,[1,3,2]);
    OMEGA = zeros(tensor_size(1),tensor_size(1),tensor_size(2));
    p = tensor_size(1);
    q = tensor_size(2);

    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(2),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end

    time = cell(q,1);
    OBJ = cell(q,1);
    parfor t = 1:q
        X_part = squeeze(X_ten_M(:,:,t));
        [ omega, ~, ~, obj,t_p] = subproblem( X_part, lambda, eta);
        OMEGA(:,:,t) = omega;
        time{t} = t_p;
        OBJ{t} = obj;
    end

    switch M
        case "fft"
            M_matrix = dftmtx(q);
        case "I"
            M_matrix = eye(q);
    end

    SCORE = zeros(p,q);
    for l = 1:p
        O_slice = squeeze(OMEGA(:,l,:));
        for h = 1:q
            m = M_matrix(:,h);
            Contri = O_slice.*repmat(m.',[p,1]);
            SCORE(l,h) = sqrt(sum(abs(Contri).^2,"all"));
        end
    end

elseif orientation == 2
    %% if direction of 2-mode is chosen
    X_ten = permute(X_ten,[1,3,2]);
    X_ten = permute(X_ten,[3,2,1]);
    OMEGA = zeros(tensor_size(2),tensor_size(2),tensor_size(1));
    p = tensor_size(2);
    q = tensor_size(1);

    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(1),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end

    time = cell(q,1);
    OBJ = cell(q,1);
    parfor t = 1:q
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~, obj,t_p] = subproblem( X_part, lambda, eta);
        OMEGA(:,:,t) = omega;
        time{t} = t_p;
        OBJ{t} = obj;
    end

    switch M
        case "fft"
            M_matrix = dftmtx(q);
        case "I"
            M_matrix = eye(q);
    end

    SCORE = zeros(p,q);
    for l = 1:p
        O_slice = squeeze(OMEGA(:,l,:));
        for h = 1:q
            m = M_matrix(:,h);
            Contri = O_slice.*repmat(m.',[p,1]);
            SCORE(l,h) = sqrt(sum(abs(Contri).^2,"all"));
        end
    end
    SCORE = SCORE.';
end

score = reshape(SCORE,nFea,1);
[~,id] = sort(score,'descend');
output.id = id;
output.score = score;
output.time_per = time;
output.OBJ = OBJ;

end

