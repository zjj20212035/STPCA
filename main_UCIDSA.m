%%%% a script to conDPct clustering experiments on feature selection
%%%% methods
clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
%% specify dataset, storage path and methods
data_name = 'UCIDSA';
data_path = strrep( pwd ,'\data','');
chosen_method = ["STPCA-MP-Dir1","STPCA-MP-Dir2"];

[save_path,keyword]= create_path(data_name,chosen_method);
load([data_path,'\data\',data_name]);
if tensor_based == 0 % 'tensor_based' is set in the loaded data and shows whether the data is a tensor or not.  
   % for non-tensor-based data, the original data matrix should be nFea*nSamp
   [nFea,nSamp] = size(data);
else 
   nSamp = length(label);
   nFea = prod(size(data,1,2));
end
%% configuration
struct.data = data;
struct.label = label;
struct.Dname = data_name;
struct.data_tensor_based = tensor_based; % whether it is tensor data, contained in the loaded dataset
struct.save_path = save_path;
struct.tensor_size = tensor_size;% the size of a tensor sample, for non-tensor data, 'tensor size' should be a [nFea,1] vector.
if tensor_based == 1
   struct.tensor_type = tensor_type; % the type of data tensor: "tube-wise",or "slice-wise"
end
grid_search = 10.^[-2:2];
NumFS = length(DiscFea); % the IDs of the features with the greatest Between-class-variance are confirmed in advance and stored in the 'DiscFea'
struct.NumFS = NumFS;
disp_method = [];
class_num = length(unique(struct.label));
for k = 1:length(chosen_method)
    if k == 1
        disp_method = [disp_method,char(chosen_method(k))];
    else
        disp_method = [disp_method,', ',char(chosen_method(k))];
    end
end

disp(['Methods: ',disp_method])

try 
    parpool;
catch
    disp('parpool already available')
end

%% Start up
%% STPCA-MP-Dir1
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1.lambda = grid_search;
    para_STPCA_MP_Dir1.eta = grid_search;
    para_special_STPCA_MP_Dir1.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1.direction = 1;
    para_special_STPCA_MP_Dir1.num_direction = 1;
    para_special_STPCA_MP_Dir1.M = "I";
    struct_STPCA_MP_Dir1 = struct;
    struct_STPCA_MP_Dir1.para = para_STPCA_MP_Dir1;
    struct_STPCA_MP_Dir1.para_special = para_special_STPCA_MP_Dir1;
    struct_STPCA_MP_Dir1.method = @STPCA_MP;
    struct_STPCA_MP_Dir1.Fname = 'STPCA-MP-Dir1';
    struct_STPCA_MP_Dir1.Ten_based = 1;
    [result_STPCA_MP_Dir1] = AlgExecution(struct_STPCA_MP_Dir1);
    save(save_path)
end

%% STPCA-MP-Dir2
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2.lambda = grid_search;
    para_STPCA_MP_Dir2.eta = grid_search;
    para_special_STPCA_MP_Dir2.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2.direction = 2;
    para_special_STPCA_MP_Dir2.num_direction = 1;
    para_special_STPCA_MP_Dir2.M = "I";
    struct_STPCA_MP_Dir2 = struct;
    struct_STPCA_MP_Dir2.para = para_STPCA_MP_Dir2;
    struct_STPCA_MP_Dir2.para_special = para_special_STPCA_MP_Dir2;
    struct_STPCA_MP_Dir2.method = @STPCA_MP;
    struct_STPCA_MP_Dir2.Fname = 'STPCA-MP-Dir2';
    struct_STPCA_MP_Dir2.Ten_based = 1;
    [result_STPCA_MP_Dir2] = AlgExecution(struct_STPCA_MP_Dir2);
    save(save_path)
end

%% Metrics computation
[mean_POC_STPCA_MP_Dir1,result_STPCA_MP_Dir1] = POC(result_STPCA_MP_Dir1,DiscFea);
[mean_POC_STPCA_MP_Dir2,result_STPCA_MP_Dir2] = POC(result_STPCA_MP_Dir2,DiscFea);


best_POC_STPCA_MP_Dir1 = result_STPCA_MP_Dir1.maximal_POC;
best_POC_STPCA_MP_Dir2 = result_STPCA_MP_Dir2.maximal_POC;

