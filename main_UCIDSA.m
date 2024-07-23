%%%% a script to conduct clustering experiments on feature selection
%%%% methods
clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
%% specify dataset, storage path and methods
data_name = 'UCIDSA';
data_path = strrep( pwd ,'\data','');
chosen_method = ["All"];

[save_path,keyword]= create_path(data_name,chosen_method);
load([data_path,'\data\',data_name]);
if tensor_based == 0 % 'tensor_based' is set in the loaded data and shows whether the data is a tensor or not.  
   data  = data'; % for non-tensor-based data, the original data matrix should be nSamp*nFea
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
   struct.tensor_type = tensor_type; % the type of data tensor: "tube-wise" or "slice-wise"
end
grid_search = 10.^[-2:2];
NumFS = length(DiscFea);
struct.NumFS = NumFS;
disp_method = [];
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
%% STPCA-DP-1SD
if  ~isempty(find(chosen_method == "STPCA-DP-1SD", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_DP_1SD.lambda = grid_search;
    para_STPCA_DP_1SD.eta = grid_search;
    para_special_STPCA_DP_1SD.S = {1};
    struct_STPCA_DP_1SD = struct;
    struct_STPCA_DP_1SD.para = para_STPCA_DP_1SD;
    struct_STPCA_DP_1SD.para_special = para_special_STPCA_DP_1SD;
    struct_STPCA_DP_1SD.method = @STPCA_DP;
    struct_STPCA_DP_1SD.Fname = 'STPCA-DP-1SD';
    struct_STPCA_DP_1SD.Ten_based = 1;
    [result_STPCA_DP_1SD] = AlgExecution(struct_STPCA_DP_1SD);
    save(save_path)
end

%% STPCA-DP-2SD
if  ~isempty(find(chosen_method == "STPCA-DP-2SD", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_DP_2SD.lambda = grid_search;
    para_STPCA_DP_2SD.eta = grid_search;
    para_special_STPCA_DP_2SD.S = {1,2};
    struct_STPCA_DP_2SD = struct;
    struct_STPCA_DP_2SD.para = para_STPCA_DP_2SD;
    struct_STPCA_DP_2SD.para_special = para_special_STPCA_DP_2SD;
    struct_STPCA_DP_2SD.method = @STPCA_DP;
    struct_STPCA_DP_2SD.Fname = 'STPCA-DP-2SD';
    struct_STPCA_DP_2SD.Ten_based = 1;
    [result_STPCA_DP_2SD] = AlgExecution(struct_STPCA_DP_2SD);
    save(save_path)
end

%% STPCA-DP-MD
if  ~isempty(find(chosen_method == "STPCA-DP-MD", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_DP_MD.lambda = grid_search;
    para_STPCA_DP_MD.eta = grid_search;
    para_special_STPCA_DP_MD.S = {[1,2]};
    struct_STPCA_DP_MD = struct;
    struct_STPCA_DP_MD.para = para_STPCA_DP_MD;
    struct_STPCA_DP_MD.para_special = para_special_STPCA_DP_MD;
    struct_STPCA_DP_MD.method = @STPCA_DP;
    struct_STPCA_DP_MD.Fname = 'STPCA-DP-MD';
    struct_STPCA_DP_MD.Ten_based = 1;
    [result_STPCA_DP_MD] = AlgExecution(struct_STPCA_DP_MD);
    save(save_path)
end

%% STPCA-MP-Dir1
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1.lambda = grid_search;
    para_STPCA_MP_Dir1.eta = grid_search;
    para_special_STPCA_MP_Dir1.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1.direction = 1;
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
    para_special_STPCA_MP_Dir2.direction = 1;
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

%% SOGFS
if  ~isempty(find(chosen_method == "SOGFS", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_SOGFS.gamma = grid_search;
    para_special_SOGFS.dim_projection = length(unique(struct.label));
    para_special_SOGFS.num_clusters = length(unique(struct.label));
    para_special_SOGFS.num_nearest_neighbors = 5;
    struct_SOGFS = struct;
    struct_SOGFS.para = para_SOGFS;
    struct_SOGFS.para_special = para_special_SOGFS;
    struct_SOGFS.method = @SOGFS;
    struct_SOGFS.Fname = 'SOGFS';
    struct_SOGFS.Ten_based = 0;
    [result_SOGFS] = AlgExecution(struct_SOGFS);
    save(save_path)
end


%% Metrics computation
POC_STPCA_DP_1SD = POC(result_STPCA_DP_1SD,DiscFea);
POC_STPCA_DP_2SD = POC(result_STPCA_DP_2SD,DiscFea);
POC_STPCA_DP_MD = POC(result_STPCA_DP_MD,DiscFea);
POC_STPCA_MP_Dir1 = POC(result_STPCA_MP_Dir1,DiscFea);
POC_STPCA_MP_Dir2 = POC(result_STPCA_MP_Dir2,DiscFea);


POTC_STPCA_DP_1SD = POTC(result_STPCA_DP_1SD,DiscFea);
POTC_STPCA_DP_2SD = POTC(result_STPCA_DP_2SD,DiscFea);
POTC_STPCA_DP_MD = POTC(result_STPCA_DP_MD,DiscFea);
POTC_STPCA_MP_Dir1 = POTC(result_STPCA_MP_Dir1,DiscFea);
POTC_STPCA_MP_Dir2 = POTC(result_STPCA_MP_Dir2,DiscFea);
save(save_path)


