clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
%% specify dataset, storage path and methods
data_name = 'pie_normalized';
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
struct.N = 30; % the repeated times of k-means 
if tensor_based == 1
   struct.tensor_type = tensor_type; % the type of data tensor: "element-wise", "tube-wise",or "slice-wise"
end
% if visualization is not one's wish, then set figpara.state = "off"
struct.figpara.state = "off";
struct.figpara.NumFS_fig = 100;
struct.figpara.type = 'jpg';
grid_search = 10.^[-2:2];
%NumFS = length(DiscFea); % the number of selected features
%NumFS = 10:20:110;
NumFS = 50:50:300;
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

%% STPCA-MP-Dir1-I
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1-I", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1_I.lambda = grid_search;
    para_STPCA_MP_Dir1_I.eta = grid_search;
    para_special_STPCA_MP_Dir1_I.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1_I.num_direction = 1;
    para_special_STPCA_MP_Dir1_I.direction = 1;
    para_special_STPCA_MP_Dir1_I.M = "I";
    struct_STPCA_MP_Dir1_I = struct;
    struct_STPCA_MP_Dir1_I.para = para_STPCA_MP_Dir1_I;
    struct_STPCA_MP_Dir1_I.para_special = para_special_STPCA_MP_Dir1_I;
    struct_STPCA_MP_Dir1_I.method = @STPCA_MP;
    struct_STPCA_MP_Dir1_I.Fname = 'STPCA-MP-Dir1-I';
    struct_STPCA_MP_Dir1_I.Ten_based = 1;
    [result_STPCA_MP_Dir1_I] = AlgExecution(struct_STPCA_MP_Dir1_I);
    [result_STPCA_MP_Dir1_I,acc_STPCA_MP_Dir1_I,nmi_STPCA_MP_Dir1_I] = ClusterExp(struct_STPCA_MP_Dir1_I,result_STPCA_MP_Dir1_I);
    save(save_path)
end

%% STPCA-MP-Dir1-fft
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1-fft", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1_fft.lambda = grid_search;
    para_STPCA_MP_Dir1_fft.eta = grid_search;
    para_special_STPCA_MP_Dir1_fft.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1_fft.num_direction = 1;
    para_special_STPCA_MP_Dir1_fft.direction = 1;
    para_special_STPCA_MP_Dir1_fft.M = "fft";
    struct_STPCA_MP_Dir1_fft = struct;
    struct_STPCA_MP_Dir1_fft.para = para_STPCA_MP_Dir1_fft;
    struct_STPCA_MP_Dir1_fft.para_special = para_special_STPCA_MP_Dir1_fft;
    struct_STPCA_MP_Dir1_fft.method = @STPCA_MP;
    struct_STPCA_MP_Dir1_fft.Fname = 'STPCA-MP-Dir1-fft';
    struct_STPCA_MP_Dir1_fft.Ten_based = 1;
    [result_STPCA_MP_Dir1_fft] = AlgExecution(struct_STPCA_MP_Dir1_fft);
    [result_STPCA_MP_Dir1_fft,acc_STPCA_MP_Dir1_fft,nmi_STPCA_MP_Dir1_fft] = ClusterExp(struct_STPCA_MP_Dir1_fft,result_STPCA_MP_Dir1_fft);
    save(save_path)
end

%% STPCA-MP-Dir1-random
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1-random", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1_random.lambda = grid_search;
    para_STPCA_MP_Dir1_random.eta = grid_search;
    para_special_STPCA_MP_Dir1_random.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1_random.num_direction = 1;
    para_special_STPCA_MP_Dir1_random.direction = 1;
    mode_2_data = tenmat(data,2);
    random_square = rand(size(mode_2_data,1),size(mode_2_data,1));
    [M_matrix,~,~] = eig(random_square*random_square');
    para_special_STPCA_MP_Dir1_random.M = "M";
    para_special_STPCA_MP_Dir1_random.M_matrix = M_matrix;
    struct_STPCA_MP_Dir1_random = struct;
    struct_STPCA_MP_Dir1_random.para = para_STPCA_MP_Dir1_random;
    struct_STPCA_MP_Dir1_random.para_special = para_special_STPCA_MP_Dir1_random;
    struct_STPCA_MP_Dir1_random.method = @STPCA_MP;
    struct_STPCA_MP_Dir1_random.Fname = 'STPCA-MP-Dir1-random';
    struct_STPCA_MP_Dir1_random.Ten_based = 1;
    [result_STPCA_MP_Dir1_random] = AlgExecution(struct_STPCA_MP_Dir1_random);
    [result_STPCA_MP_Dir1_random,acc_STPCA_MP_Dir1_random,nmi_STPCA_MP_Dir1_random] = ClusterExp(struct_STPCA_MP_Dir1_random,result_STPCA_MP_Dir1_random);
    save(save_path)
end

%% STPCA-MP-Dir1-mode
if  ~isempty(find(chosen_method == "STPCA-MP-Dir1-mode", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir1_mode.lambda = grid_search;
    para_STPCA_MP_Dir1_mode.eta = grid_search;
    para_special_STPCA_MP_Dir1_mode.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir1_mode.num_direction = 1;
    para_special_STPCA_MP_Dir1_mode.direction = 1;
    mode_2_data = tenmat(data,2);
    [M_matrix,S,~] = eig(double(mode_2_data*mode_2_data'));
    para_special_STPCA_MP_Dir1_mode.M = "M";
    para_special_STPCA_MP_Dir1_mode.M_matrix = M_matrix;
    struct_STPCA_MP_Dir1_mode = struct;
    struct_STPCA_MP_Dir1_mode.para = para_STPCA_MP_Dir1_mode;
    struct_STPCA_MP_Dir1_mode.para_special = para_special_STPCA_MP_Dir1_mode;
    struct_STPCA_MP_Dir1_mode.method = @STPCA_MP;
    struct_STPCA_MP_Dir1_mode.Fname = 'STPCA-MP-Dir1-mode';
    struct_STPCA_MP_Dir1_mode.Ten_based = 1;
    [result_STPCA_MP_Dir1_mode] = AlgExecution(struct_STPCA_MP_Dir1_mode);
    [result_STPCA_MP_Dir1_mode,acc_STPCA_MP_Dir1_mode,nmi_STPCA_MP_Dir1_mode] = ClusterExp(struct_STPCA_MP_Dir1_mode,result_STPCA_MP_Dir1_mode);
    save(save_path)
end

%% STPCA-MP-Dir2-I
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2-I", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2_I.lambda = grid_search;
    para_STPCA_MP_Dir2_I.eta = grid_search;
    para_special_STPCA_MP_Dir2_I.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2_I.num_direction = 1;
    para_special_STPCA_MP_Dir2_I.direction = 2;
    para_special_STPCA_MP_Dir2_I.M = "I";
    struct_STPCA_MP_Dir2_I = struct;
    struct_STPCA_MP_Dir2_I.para = para_STPCA_MP_Dir2_I;
    struct_STPCA_MP_Dir2_I.para_special = para_special_STPCA_MP_Dir2_I;
    struct_STPCA_MP_Dir2_I.method = @STPCA_MP;
    struct_STPCA_MP_Dir2_I.Fname = 'STPCA-MP-Dir2-I';
    struct_STPCA_MP_Dir2_I.Ten_based = 1;
    [result_STPCA_MP_Dir2_I] = AlgExecution(struct_STPCA_MP_Dir2_I);
    [result_STPCA_MP_Dir2_I,acc_STPCA_MP_Dir2_I,nmi_STPCA_MP_Dir2_I] = ClusterExp(struct_STPCA_MP_Dir2_I,result_STPCA_MP_Dir2_I);  
    save(save_path)
end

%% STPCA-MP-Dir2-fft
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2-fft", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2_fft.lambda = grid_search;
    para_STPCA_MP_Dir2_fft.eta = grid_search;
    para_special_STPCA_MP_Dir2_fft.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2_fft.num_direction = 1;
    para_special_STPCA_MP_Dir2_fft.direction = 2;
    para_special_STPCA_MP_Dir2_fft.M = "fft";
    struct_STPCA_MP_Dir2_fft = struct;
    struct_STPCA_MP_Dir2_fft.para = para_STPCA_MP_Dir2_fft;
    struct_STPCA_MP_Dir2_fft.para_special = para_special_STPCA_MP_Dir2_fft;
    struct_STPCA_MP_Dir2_fft.method = @STPCA_MP;
    struct_STPCA_MP_Dir2_fft.Fname = 'STPCA-MP-Dir2-fft';
    struct_STPCA_MP_Dir2_fft.Ten_based = 1;
    [result_STPCA_MP_Dir2_fft] = AlgExecution(struct_STPCA_MP_Dir2_fft);
    [result_STPCA_MP_Dir2_fft,acc_STPCA_MP_Dir2_fft,nmi_STPCA_MP_Dir2_fft] = ClusterExp(struct_STPCA_MP_Dir2_fft,result_STPCA_MP_Dir2_fft);
    save(save_path)
end

%% STPCA-MP-Dir2-random
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2-random", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2_random.lambda = grid_search;
    para_STPCA_MP_Dir2_random.eta = grid_search;
    para_special_STPCA_MP_Dir2_random.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2_random.num_direction = 1;
    para_special_STPCA_MP_Dir2_random.direction = 2;
    mode_1_data = tenmat(data,1);
    random_square = rand(size(mode_1_data,1),size(mode_1_data,1));
    [M_matrix,~,~] = eig(random_square*random_square');
    para_special_STPCA_MP_Dir2_random.M = "M";
    para_special_STPCA_MP_Dir2_random.M_matrix = M_matrix;
    struct_STPCA_MP_Dir2_random = struct;
    struct_STPCA_MP_Dir2_random.para = para_STPCA_MP_Dir2_random;
    struct_STPCA_MP_Dir2_random.para_special = para_special_STPCA_MP_Dir2_random;
    struct_STPCA_MP_Dir2_random.method = @STPCA_MP;
    struct_STPCA_MP_Dir2_random.Fname = 'STPCA-MP-Dir2-random';
    struct_STPCA_MP_Dir2_random.Ten_based = 1;
    [result_STPCA_MP_Dir2_random] = AlgExecution(struct_STPCA_MP_Dir2_random);
    [result_STPCA_MP_Dir2_random,acc_STPCA_MP_Dir2_random,nmi_STPCA_MP_Dir2_random] = ClusterExp(struct_STPCA_MP_Dir2_random,result_STPCA_MP_Dir2_random); 
    save(save_path)
end

%% STPCA-MP-Dir2-mode
if  ~isempty(find(chosen_method == "STPCA-MP-Dir2-mode", 1)) || ~isempty(find(chosen_method == "All", 1))
    para_STPCA_MP_Dir2_mode.lambda = grid_search;
    para_STPCA_MP_Dir2_mode.eta = grid_search;
    para_special_STPCA_MP_Dir2_mode.imag_size = struct.tensor_size;
    para_special_STPCA_MP_Dir2_mode.num_direction = 1;
    para_special_STPCA_MP_Dir2_mode.direction = 2;
    mode_1_data = tenmat(data,1);
    [M_matrix,~,~] = eig(double(mode_1_data*mode_1_data'));
    para_special_STPCA_MP_Dir2_mode.M = "M";
    para_special_STPCA_MP_Dir2_mode.M_matrix = M_matrix;
    struct_STPCA_MP_Dir2_mode = struct;
    struct_STPCA_MP_Dir2_mode.para = para_STPCA_MP_Dir2_mode;
    struct_STPCA_MP_Dir2_mode.para_special = para_special_STPCA_MP_Dir2_mode;
    struct_STPCA_MP_Dir2_mode.method = @STPCA_MP;
    struct_STPCA_MP_Dir2_mode.Fname = 'STPCA-MP-Dir2-mode';
    struct_STPCA_MP_Dir2_mode.Ten_based = 1;
    [result_STPCA_MP_Dir2_mode] = AlgExecution(struct_STPCA_MP_Dir2_mode);
    [result_STPCA_MP_Dir2_mode,acc_STPCA_MP_Dir2_mode,nmi_STPCA_MP_Dir2_mode] = ClusterExp(struct_STPCA_MP_Dir2_mode,result_STPCA_MP_Dir2_mode);  
    save(save_path)
end



