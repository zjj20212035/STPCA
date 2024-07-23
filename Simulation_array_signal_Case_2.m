%%%% Simulation
clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
addpath(genpath([pwd,'\result\analysis']))
%% specify dataset, storage path and methods
keyword = 'Vertical';
data_path = [pwd,'\data\Array signal\',keyword,'\'];
files = dir(fullfile(data_path,'*.mat'));
num_mat = 0;
name = [];
for k = 1:length(files)
    num_mat = num_mat + 1;
    name = [name,string(files(k).name)];
end

try 
    parpool;
catch
    disp('parpool already available')
end

POC_STPCA_DP_2SD = zeros(num_mat,1);
POC_STPCA_DP_MD = zeros(num_mat,1);
POC_STPCA_MP_Dir1 = zeros(num_mat,1);
POC_STPCA_MP_Dir2 = zeros(num_mat,1);

POTC_STPCA_DP_2SD = zeros(num_mat,1);
POTC_STPCA_DP_MD = zeros(num_mat,1);
POTC_STPCA_MP_Dir1 = zeros(num_mat,1);
POTC_STPCA_MP_Dir2 = zeros(num_mat,1);

for p = 1:num_mat
    data_name = [keyword,'_',num2str(p)];
    chosen_method = ["All"];
    [save_path,~]= create_path(data_name,chosen_method);
    load([data_path,char(name(p))]);
    if tensor_based == 0 % 'tensor_based' is set in the loaded data and shows whether the data is a tensor or not.
        data = data';
        [nFea,nSamp] = size(data);
    else
        nSamp = length(label);
        nFea = prod(size(data,1,2));
    end
    %% configuration
    grid_search = 10.^[-2:2];
    struct.data = data;
    struct.label = label;
    struct.Dname = data_name;
    struct.data_tensor_based = tensor_based; % whether it is tensor data, containing in the loaded dataset
    struct.save_path = save_path;
    struct.tensor_size = tensor_size;% the type of data tensor: "tube-wise" or "slice-wise"
    if tensor_based == 1
        struct.tensor_type = tensor_type;
    end
    NumFS = length(error_position);
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

    %% Start up
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
        para_special_STPCA_MP_Dir2.direction = 2;
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
    %% Metric analysis
    POC_STPCA_DP_2SD(p) = POC(result_STPCA_DP_2SD,error_position);
    POC_STPCA_DP_MD(p) = POC(result_STPCA_DP_MD,error_position);
    POC_STPCA_MP_Dir1(p) = POC(result_STPCA_MP_Dir1,error_position);
    POC_STPCA_MP_Dir2(p) = POC(result_STPCA_MP_Dir2,error_position);

    POTC_STPCA_DP_2SD(p) = POTC(result_STPCA_DP_2SD,error_position);
    POTC_STPCA_DP_MD(p) = POTC(result_STPCA_DP_MD,error_position);
    POTC_STPCA_MP_Dir1(p) = POTC(result_STPCA_MP_Dir1,error_position);
    POTC_STPCA_MP_Dir2(p) = POTC(result_STPCA_MP_Dir2,error_position);
    save(save_path)
end

[mean_POC_STPCA_DP_2SD, std_POC_STPCA_DP_2SD] = mean_std(POC_STPCA_DP_2SD);
[mean_POC_STPCA_DP_MD, std_POC_STPCA_DP_MD] = mean_std(POC_STPCA_DP_MD);
[mean_POC_STPCA_MP_Dir1, std_POC_STPCA_MP_Dir1] = mean_std(POC_STPCA_MP_Dir1);
[mean_POC_STPCA_MP_Dir2, std_POC_STPCA_MP_Dir2] = mean_std(POC_STPCA_MP_Dir2);

[mean_POTC_STPCA_DP_2SD, std_POTC_STPCA_DP_2SD] = mean_std(POTC_STPCA_DP_2SD);
[mean_POTC_STPCA_DP_MD, std_POTC_STPCA_DP_MD] = mean_std(POTC_STPCA_DP_MD);
[mean_POTC_STPCA_MP_Dir1, std_POTC_STPCA_MP_Dir1] = mean_std(POTC_STPCA_MP_Dir1);
[mean_POTC_STPCA_MP_Dir2, std_POTC_STPCA_MP_Dir2] = mean_std(POTC_STPCA_MP_Dir2);

