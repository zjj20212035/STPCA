function [output] = STPCA_MP(X_ten,para,para_special)
%Sparse Tensor PCA based on *M product
%X_ten: a data tensor of size fea1*fea2*Sample
% make sure parpool is opened in advance!!
lambda = para.lambda;
eta = para.eta;
imag_size = para_special.imag_size;
M = para_special.M;

nFea = prod(imag_size);
size_X = size(X_ten);
direction = para_special.direction;

if direction == 1
    %% if direction 1 is chosen
    X_ten = permute(X_ten,[1,3,2]);
    OMEGA = zeros(imag_size(1),imag_size(1),imag_size(2));
    % *M product
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(2),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M,3));
    end
    % slice-by-slice 
    parfor t = 1:imag_size(2)
        X_part = squeeze(X_ten_M(:,:,t));
        [ omega, ~, ~, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA(:,:,t) = omega;
    end
    % *M product
    switch M
        case "fft"
            OMEGA = ifft(OMEGA, imag_size(2), 3);
        case "I"
            OMEGA = OMEGA;
        otherwise
            OMEGA = double(ttm(tensor(OMEGA),M^-1,3));
    end
    SCORE = squeeze(sum(abs(OMEGA).^2));

elseif direction == 2
    %% if direction 2 is chosen
    X_ten = permute(X_ten,[1,3,2]);
    X_ten = permute(X_ten,[3,2,1]);
    OMEGA = zeros(imag_size(2),imag_size(2),imag_size(1));
    % *M product
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(1),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M,3));
    end
    % slice-by-slice
    parfor t = 1:imag_size(1)
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~ ] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA(:,:,t) = omega;
    end
    % *M product
    switch M
        case "fft"
            OMEGA = ifft(OMEGA, imag_size(1), 3);
        case "I"
            OMEGA = OMEGA;
        otherwise
            OMEGA = double(ttm(tensor(OMEGA),M^-1,3));
    end
    SCORE = squeeze(sum(abs(OMEGA).^2)).';
end


score = reshape(SCORE,nFea,1);
[~,id] = sort(score,'descend');
output.id = id;
output.score = score;

end

