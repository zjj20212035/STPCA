function [output] = STPCA_MP(X_ten,para,para_special)
%Sparse Tensor PCA based on *M product
%X_ten: a data tensor that contains a series of 2D samples
lambda = para.lambda;
eta = para.eta;
imag_size = para_special.imag_size;
num_direction = para_special.num_direction;
M = para_special.M;
if M~="fft" && M~="I"
   M_matrix = para_special.M_matrix;
end

nFea = prod(imag_size);
size_X = size(X_ten);
nSamp = size_X(end);
num_diag = sum(imag_size)-3;
first_diag = -(imag_size(1)-2);
max_length_diag = min(imag_size);
diag_index =  first_diag:(first_diag+num_diag-1);
%diag_para_shrink_ratio = 0.9;

if num_direction == 1
    direction = para_special.direction;
    if direction == 1
        %% if direction 1 is chosen （D^o = {1,3,2}）
        X_ten = permute(X_ten,[1,3,2]);
        OMEGA = zeros(imag_size(1),imag_size(1),imag_size(2));
        p = imag_size(1);
        q = imag_size(2);
        switch M
            case "fft"
                X_ten_M = fft(X_ten,size_X(2),3);
            case "I"
                X_ten_M = X_ten;
            otherwise
                X_ten_M = tensor(X_ten);
                X_ten_M = double(ttm(X_ten_M,M_matrix,3));
        end

        parfor t = 1:q
            X_part = squeeze(X_ten_M(:,:,t));
            [ omega, ~, ~, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
            OMEGA(:,:,t) = omega;
        end

        
        switch M
            case "fft"
                M_matrix = dftmtx(q);
            case "I"
                M_matrix = eye(q);
            otherwise
                M_matrix = M_matrix;
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
        
    elseif direction == 2
        %% if direction 2 is chosen （D^o = {2,3,1}）
        X_ten = permute(X_ten,[1,3,2]);
        X_ten = permute(X_ten,[3,2,1]);
        OMEGA = zeros(imag_size(2),imag_size(2),imag_size(1));
        p = imag_size(2);
        q = imag_size(1);        
        switch M
            case "fft"
                X_ten_M = fft(X_ten,size_X(1),3);
            case "I"
                X_ten_M = X_ten;
            otherwise
                X_ten_M = tensor(X_ten);
                X_ten_M = double(ttm(X_ten_M,M_matrix,3));
        end
        parfor t = 1:q
            X_part = squeeze(X_ten_M(:,:,t));
            [omega, ~, ~ ] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
            OMEGA(:,:,t) = omega;
        end

        switch M
            case "fft"
                M_matrix = dftmtx(q);
            case "I"
                M_matrix = eye(q);
            otherwise
                M_matrix = M_matrix;
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

elseif num_direction == 2 % simultaneouly considering two directions (still waiting for further research)
    OMEGA_1 = zeros(imag_size(1),imag_size(1),imag_size(2));
    OMEGA_2 = zeros(imag_size(2),imag_size(2),imag_size(1));
    %% direction 1
    X_ten = permute(X_ten,[1,3,2]);
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(2),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end
    parfor t = 1:imag_size(2)
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~, ~ ] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA_1(:,:,t) = omega;
    end
    switch M
        case "fft"
            OMEGA_1 = ifft(OMEGA_1, imag_size(2), 3);
        case "I"
            OMEGA_1 = OMEGA_1;
        otherwise
            OMEGA_1 = double(ttm(tensor(OMEGA_1),M_matrix^-1,3));
    end
    score_1 = squeeze(sum(abs(OMEGA_1).^2));
   
    %% direction 2
    X_ten = permute(X_ten,[3,2,1]);
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(1),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end
    parfor t = 1:imag_size(1)
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~ ] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA_2(:,:,t) = omega;
    end
    switch M
        case "fft"
            OMEGA_2 = ifft(OMEGA_2, imag_size(1), 3);
        case "I"
            OMEGA_2 = OMEGA_2;
        otherwise
            OMEGA_2 = double(ttm(tensor(OMEGA_2),M_matrix^-1,3));
    end
    score_2 = squeeze(sum(abs(OMEGA_2).^2));
    SCORE = score_1 + score_2.' ;


elseif num_direction == 4
    score_3 = zeros(imag_size);
    score_4 = zeros(imag_size);
    OMEGA_1 = zeros(imag_size(1),imag_size(1),imag_size(2));
    OMEGA_2 = zeros(imag_size(2),imag_size(2),imag_size(1));
    Bin_1 = zeros(max_length_diag,num_diag); % to store the elements of each dianonal slice
    Bin_2 = zeros(max_length_diag,num_diag);
    flag = 1;
    %% direction 1
    X_ten = permute(X_ten,[1,3,2]);
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(2),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end
    parfor t = 1:imag_size(2)
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA_1(:,:,t) = omega;
    end
    switch M
        case "fft"
            OMEGA_1 = ifft(OMEGA_1, imag_size(2), 3);
        case "I"
            OMEGA_1 = OMEGA_1;
        otherwise
            OMEGA_1 = double(ttm(tensor(OMEGA_1),M_matrix^-1,3));
    end
    score_1 = squeeze(sum(abs(OMEGA_1).^2));
   
    %% direction 2
    X_ten = permute(X_ten,[3,2,1]);
    switch M
        case "fft"
            X_ten_M = fft(X_ten,size_X(1),3);
        case "I"
            X_ten_M = X_ten;
        otherwise
            X_ten_M = tensor(X_ten);
            X_ten_M = double(ttm(X_ten_M,M_matrix,3));
    end
    parfor t = 1:imag_size(1)
        X_part = squeeze(X_ten_M(:,:,t));
        [omega, ~, ~, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        OMEGA_2(:,:,t) = omega;
    end
    switch M
        case "fft"
            OMEGA_2 = ifft(OMEGA_2, imag_size(1), 3);
        case "I"
            OMEGA_2 = OMEGA_2;
        otherwise
            OMEGA_2 = double(ttm(tensor(OMEGA_2),M_matrix^-1,3));
    end
    score_2 = squeeze(sum(abs(OMEGA_2).^2));

    %% direction 3
    X_ten = permute(X_ten,[3,2,1]);
    X_ten = permute(X_ten,[1,3,2]);
    for t = diag_index
        length_diag = length(diag(X_ten(:,:,1),t));
        X_part = zeros(length_diag,nSamp);
        if t <= 0
            row_index = (1-t):imag_size(1);
            column_index = 1:length(row_index);
        else
            column_index = (1+t):imag_size(2);
            row_index = 1:length(column_index);
        end
        for n = 1:length_diag
            X_part(n,:) = squeeze(X_ten(row_index(n),column_index(n),:));
        end
        [~, ~, score, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        if imag_size(1) >= imag_size(2)
            if t<=0
                Bin_1(1:length_diag,flag) = score;
            else
                Bin_1((end-length_diag+1):end, flag) = score;
            end
        else
            if t<=0
                Bin_1( (end-length_diag+1):end, flag) = score;
            else
                Bin_1(1:length_diag, flag) = score;
            end
        end
        flag = flag + 1;
    end
    score_3 = full(spdiags(Bin_1,diag_index,score_3));


    %% direction 4
    X_ten_fliplr = fliplr(X_ten);
    flag = 1;
    for t = diag_index
        length_diag = length(diag(X_ten(:,:,1),t));
        X_part = zeros(length_diag,nSamp);
        if t <= 0
            row_index = (1-t):imag_size(1);
            column_index = 1:length(row_index);
        else
            column_index = (1+t):imag_size(2);
            row_index = 1:length(column_index);
        end
        for n = 1:length_diag
            X_part(n,:) = squeeze(X_ten_fliplr(row_index(n),column_index(n),:));
        end
        [~, ~, score, ~] = SPCA_PSD_for_sequential_slice( X_part, lambda, eta);
        if imag_size(1) >= imag_size(2)
            if t<=0
                Bin_2(1:length_diag,flag) = score;
            else
                Bin_2( (end-length_diag+1):end, flag) = score;
            end
        else
            if t<=0
                Bin_2( (end-length_diag+1):end, flag) = score;
            else
                Bin_2(1:length_diag,flag) = score;
            end
        end
        flag = flag + 1;
    end
    score_4 = full(spdiags(Bin_2,diag_index,score_4));
    score_4 = fliplr(score_4);
    SCORE = score_1 + score_2.' + score_3 + score_4;
end


score = reshape(SCORE,nFea,1);
[~,id] = sort(score,'descend');
output.id = id;
output.score = score;

end

