%%%%%%%% configuration %%%%%%%%%
chosen_method = ["STPCA-MP"];
% create save path
[save_path,keyword]= create_path(data_name,chosen_method);
load(data_path);
if tensor_based == 0 % 0:non-tensor data;  1:tensor data
   % for non-tensor-based data, the original data matrix should be nFea*nSamp
   [nFea,nSamp] = size(data);
else 
   nSamp = length(label);
   nFea = prod(size(data,1,2));
end
%% data preprocessing
data = data/max(abs(data),[],'all'); % normalized to [-1,1]
%% experiment settings for all datasets
struct.data = data;
struct.label = label;
struct.Dname = data_name;
struct.data_tensor_based = tensor_based; % whether it is tensor data, contained in the loaded dataset
struct.save_path = save_path;
struct.tensor_size = tensor_size;% the size of a tensor sample, for non-tensor data, 'tensor size' should be a [nFea,1] vector.
struct.N = 30; % the repeated times of k-means

if tensor_based == 1
   struct.tensor_type = tensor_type; % the type of data tensor: "tube-wise" or "slice-wise"
end

%% grid-search strategy
grid_search = 10.^[-2:2]; % on all datasets

%% settings of baseline: All Features 
struct_All_Features = struct; 
struct_All_Features.NumFS = nFea;
para_All_Features.nFea = nFea;
para_All_Features.nSamp = nSamp;
struct_All_Features.para = para_All_Features;
struct_All_Features.para_special = [];
struct_All_Features.method = @All_Features;
struct_All_Features.Fname = 'All-Features';
struct_All_Features.Ten_based = 0;

%% settings of STPCA-MP
para_STPCA_MP.lambda = grid_search;
para_STPCA_MP.eta = grid_search;
para_special_STPCA_MP.tensor_size = struct.tensor_size;
para_special_STPCA_MP.M = "I";
struct_STPCA_MP = struct;
struct_STPCA_MP.para = para_STPCA_MP;
struct_STPCA_MP.method = @STPCA_MP;
struct_STPCA_MP.Fname = 'STPCA-MP';
struct_STPCA_MP.Ten_based = 1;
