%%%% a script to conDPct clustering experiments on feature selection
%%%% methods
clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
%% specify dataset, storage path and methods
data_name = 'pie_normalized';
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
struct.N = 30; % the repeated times of k-means 
if tensor_based == 1
   struct.tensor_type = tensor_type; % the type of data tensor: "element-wise", "tube-wise",or "slice-wise"
end
% if visualization is not one's wish, then set figpara.state = "off"
struct.figpara.state = "off";
struct.figpara.NumFS_fig = 100;
struct.figpara.type = 'jpg';
grid_search = 10.^[-2:2];
NumFS = 50:50:300; % the number of selected features
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

disp(['Methods: ',disp_method,';   Visualization: ', char(struct.figpara.state)])

try 
    parpool;
catch
    disp('parpool already available')
end


%% Start up
%% All features
if  ~isempty(find(chosen_method == "All features", 1))  || ~isempty(find(chosen_method == "All", 1))
    N = struct.N;
    label = struct.label;
    ACC_All_Feature = cell(1,N);
    NMI_All_Feature = cell(1,N);
    if tensor_based == 0
        X_r = data';
    else
        X_r = reshape(data,[prod(tensor_size),nSamp]);
    end
    if ~isreal(data)
        X_r = [abs(X_r);angle(X_r)];
    end
    parfor i1 = 1:N
        idx = kmeans(X_r',class_num);
        final_idx = BestMapping(label,idx );
        ACC_All_Feature{i1} = 1 - length(find(final_idx -label)) / nSamp;
        NMI_All_Feature{i1} = nmi(label, final_idx);
    end
    ACC_All_Feature = repmat(mean(cell2mat(ACC_All_Feature)),[1,length(NumFS)]);
    NMI_All_Feature = repmat(mean(cell2mat(NMI_All_Feature)),[1,length(NumFS)]);
    STD_ACC_All_Feature = std(ACC_All_Feature);
    STD_NMI_All_Feature = std(NMI_All_Feature);
    save(save_path);
end

%% STPCA-MP-Dir1
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1.lambda = grid_search;
    para_STPCA_MP_Dir1.eta = grid_search;
    para_special_STPCA_MP_Dir1.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1.num_direction = 1;
    para_special_STPCA_MP_Dir1.direction = 1;
    para_special_STPCA_MP_Dir1.M = "I";
    struct_STPCA_MP_Dir1 = struct;
    struct_STPCA_MP_Dir1.para = para_STPCA_MP_Dir1;
    struct_STPCA_MP_Dir1.para_special = para_special_STPCA_MP_Dir1;
    struct_STPCA_MP_Dir1.method = @STPCA_MP;
    struct_STPCA_MP_Dir1.Fname = 'STPCA-MP-Dir1';
    struct_STPCA_MP_Dir1.Ten_based = 1;
    [result_STPCA_MP_Dir1] = AlgExecution(struct_STPCA_MP_Dir1);
    [result_STPCA_MP_Dir1,acc_STPCA_MP_Dir1,nmi_STPCA_MP_Dir1] = ClusterExp(struct_STPCA_MP_Dir1,result_STPCA_MP_Dir1);
    save(save_path)
end

%% STPCA-MP-Dir2
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2.lambda = grid_search;
    para_STPCA_MP_Dir2.eta = grid_search;
    para_special_STPCA_MP_Dir2.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2.num_direction = 1;
    para_special_STPCA_MP_Dir2.direction = 2;
    para_special_STPCA_MP_Dir2.M = "I";
    struct_STPCA_MP_Dir2 = struct;
    struct_STPCA_MP_Dir2.para = para_STPCA_MP_Dir2;
    struct_STPCA_MP_Dir2.para_special = para_special_STPCA_MP_Dir2;
    struct_STPCA_MP_Dir2.method = @STPCA_MP;
    struct_STPCA_MP_Dir2.Fname = 'STPCA-MP-Dir2';
    struct_STPCA_MP_Dir2.Ten_based = 1;
    [result_STPCA_MP_Dir2] = AlgExecution(struct_STPCA_MP_Dir2);
    [result_STPCA_MP_Dir2,acc_STPCA_MP_Dir2,nmi_STPCA_MP_Dir2] = ClusterExp(struct_STPCA_MP_Dir2,result_STPCA_MP_Dir2);
    save(save_path)
end

