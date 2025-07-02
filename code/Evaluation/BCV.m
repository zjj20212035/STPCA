function [BCV, DiscFea] = BCV(data,label,tensor_type, NumFS)
% compute Between-class variance (BCV) of tensor data
% NumFS: number of selected feature
% DiscFea: discriminatve features
data_centralized = data - mean(data,3);
num_class = length(unique(label));
mean_class = zeros(size(data,1),size(data,2),num_class);

if tensor_type == "slice-wise"
    for h = unique(min(label)):num_class
        [a,~] = find(label==h);
        mean_class (:,:,h) = mean(data_centralized(:,:,a),3);
    end
    BCV = sum(var(mean_class,0,3),2);
elseif tensor_type == "tube-wise" || tensor_type == "element-wise"
    for h = unique(min(label)):num_class
        [a,~] = find(label==h);
        mean_class (:,:,h) = mean(data_centralized(:,:,a),3);
    end
   BCV = reshape(var(mean_class,0,3),[prod(tensor_size),1]);
end

[~,max_BCV_feature] = sort(BCV,'descend');
DiscFea = max_BCV_feature(1:NumFS);
end

