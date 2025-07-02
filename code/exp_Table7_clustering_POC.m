%%%% A reproduction of Table 7 (clustering ACC, NMI, and SS; POC) in 
%%%% Orientation-Aware Sparse Tensor PCA for Efficient Unsupervised Feature
%%%% Selection%

%% preparing environment
clear
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
addpath(genpath([pwd,'\utils']))
% start parallel computation toolbox
try 
    parpool;
catch
    disp('parpool already available')
end

%% experiment

%% tube-wise datasets
%%%%%%%%%%%% PIE %%%%%%%%%%%%%%%%%
data_name = 'PIE';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_PIE] = AlgExecution(struct_All_Features);
[result_All_Features_PIE,acc_All_Features_pie,nmi_All_Features_pie,ss_All_Features_pie,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_PIE);
% conduct STPCA-MP
para_special_STPCA_MP_PIE = para_special_STPCA_MP;
para_special_STPCA_MP_PIE.orientation = 1;
struct_STPCA_MP_PIE = struct_STPCA_MP;
struct_STPCA_MP_PIE.NumFS = 50:50:300;
struct_STPCA_MP_PIE.para_special = para_special_STPCA_MP_PIE;
[result_STPCA_MP_PIE] = AlgExecution(struct_STPCA_MP_PIE);
[result_STPCA_MP_PIE,acc_STPCA_MP_PIE,nmi_STPCA_MP_PIE,ss_STPCA_MP_PIE,...,
    ] = ClusterExp(struct_STPCA_MP_PIE,result_STPCA_MP_PIE);
save(save_path)

%%%%%%%%%%%% COIL20 %%%%%%%%%%%%%%%%%
data_name = 'COIL20';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_COIL20] = AlgExecution(struct_All_Features);
[result_All_Features_COIL20,acc_All_Features,nmi_All_Features_COIL20,ss_All_Features_COIL20,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_COIL20);
% conduct STPCA-MP
para_special_STPCA_MP_COIL20 = para_special_STPCA_MP;
para_special_STPCA_MP_COIL20.orientation = 1;
struct_STPCA_MP_COIL20 = struct_STPCA_MP;
struct_STPCA_MP_COIL20.NumFS = 50:50:300;
struct_STPCA_MP_COIL20.para_special = para_special_STPCA_MP_COIL20;
[result_STPCA_MP_COIL20] = AlgExecution(struct_STPCA_MP_COIL20);
[result_STPCA_MP_COIL20,acc_STPCA_MP_COIL20,nmi_STPCA_MP_COIL20,ss_STPCA_MP_COIL20,...,
    ] = ClusterExp(struct_STPCA_MP_COIL20,result_STPCA_MP_COIL20);
save(save_path)

%%%%%%%%%%%% JAFFE %%%%%%%%%%%%%%%%%
data_name = 'JAFFE';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_JAFFE] = AlgExecution(struct_All_Features);
[result_All_Features_JAFFE,acc_All_Features_JAFFE,nmi_All_Features_JAFFE,ss_All_Features_JAFFE,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_JAFFE);
% conduct STPCA-MP
para_special_STPCA_MP_JAFFE = para_special_STPCA_MP;
para_special_STPCA_MP_JAFFE.orientation = 1;
struct_STPCA_MP_JAFFE = struct_STPCA_MP;
struct_STPCA_MP_JAFFE.NumFS = 50:50:300;
struct_STPCA_MP_JAFFE.para_special = para_special_STPCA_MP_JAFFE;
[result_STPCA_MP_JAFFE] = AlgExecution(struct_STPCA_MP_JAFFE);
[result_STPCA_MP_JAFFE,acc_STPCA_MP_JAFFE,nmi_STPCA_MP_JAFFE,ss_STPCA_MP_JAFFE,...,
    ] = ClusterExp(struct_STPCA_MP_JAFFE,result_STPCA_MP_JAFFE);
save(save_path)

%%%%%%%%%%%% Imm40 %%%%%%%%%%%%%%%%%
data_name = 'Imm40';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_Imm40] = AlgExecution(struct_All_Features);
[result_All_Features_Imm40,acc_All_Features_Imm40,nmi_All_Features_Imm40,ss_All_Features_Imm40,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_Imm40);
% conduct STPCA-MP
para_special_STPCA_MP_Imm40 = para_special_STPCA_MP;
para_special_STPCA_MP_Imm40.orientation = 2;
struct_STPCA_MP_Imm40 = struct_STPCA_MP;
struct_STPCA_MP_Imm40.NumFS = 50:50:300;
struct_STPCA_MP_Imm40.para_special = para_special_STPCA_MP_Imm40;
[result_STPCA_MP_Imm40] = AlgExecution(struct_STPCA_MP_Imm40);
[result_STPCA_MP_Imm40,acc_STPCA_MP_Imm40,nmi_STPCA_MP_Imm40,ss_STPCA_MP_Imm40,...,
    ] = ClusterExp(struct_STPCA_MP_Imm40,result_STPCA_MP_Imm40);
save(save_path)

%%%%%%%%%%%% BreastMNIST %%%%%%%%%%%%%%%%%
data_name = 'BreastMNIST';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_BreastMNIST] = AlgExecution(struct_All_Features);
[result_All_Features_BreastMNIST,acc_All_Features_BreastMNIST,nmi_All_Features_BreastMNIST,ss_All_Features_BreastMNIST,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_BreastMNIST);
% conduct STPCA-MP
para_special_STPCA_MP_BreastMNIST = para_special_STPCA_MP;
para_special_STPCA_MP_BreastMNIST.orientation = 1;
struct_STPCA_MP_BreastMNIST = struct_STPCA_MP;
struct_STPCA_MP_BreastMNIST.NumFS = 50:50:300;
struct_STPCA_MP_BreastMNIST.para_special = para_special_STPCA_MP_BreastMNIST;
[result_STPCA_MP_BreastMNIST] = AlgExecution(struct_STPCA_MP_BreastMNIST);
[result_STPCA_MP_BreastMNIST,acc_STPCA_MP_BreastMNIST,nmi_STPCA_MP_BreastMNIST,ss_STPCA_MP_BreastMNIST,...,
    ] = ClusterExp(struct_STPCA_MP_BreastMNIST,result_STPCA_MP_BreastMNIST);
save(save_path)

%%%%%%%%%%%% USPS %%%%%%%%%%%%%%%%%
data_name = 'USPS';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% conduct All Features as baseline
[result_All_Features_USPS] = AlgExecution(struct_All_Features);
[result_All_Features_USPS,acc_All_Features_USPS,nmi_All_Features_USPS,ss_All_Features_USPS,...,
    ] = ClusterExp(struct_All_Features,result_All_Features_USPS);
% conduct STPCA-MP
para_special_STPCA_MP_USPS = para_special_STPCA_MP;
para_special_STPCA_MP_USPS.orientation = 2;
struct_STPCA_MP_USPS = struct_STPCA_MP;
struct_STPCA_MP_USPS.NumFS = 10:20:110;
struct_STPCA_MP_USPS.para_special = para_special_STPCA_MP_USPS;
[result_STPCA_MP_USPS] = AlgExecution(struct_STPCA_MP_USPS);
[result_STPCA_MP_USPS,acc_STPCA_MP_USPS,nmi_STPCA_MP_USPS,ss_STPCA_MP_USPS,...,
    ] = ClusterExp(struct_STPCA_MP_USPS,result_STPCA_MP_USPS);
save(save_path)


%% slice-wise datasets
%%%%%%%%%%%% UCIDSA %%%%%%%%%%%%%%%%%
data_name = 'UCIDSA';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% obtain features with the greatest Between-class variance (BCV) of data
NumFS = 15;
[~,DiscFea] = BCV(data, label, tensor_type, NumFS);
% conduct STPCA-MP
para_special_STPCA_MP_UCIDSA = para_special_STPCA_MP;
para_special_STPCA_MP_UCIDSA.orientation = 1;
struct_STPCA_MP_UCIDSA = struct_STPCA_MP;
struct_STPCA_MP_UCIDSA.NumFS = NumFS;
struct_STPCA_MP_UCIDSA.para_special = para_special_STPCA_MP_UCIDSA;
[result_STPCA_MP_UCIDSA] = AlgExecution(struct_STPCA_MP_UCIDSA);
[mean_POC_STPCA_MP_UCIDSA,result_STPCA_MP_UCIDSA] = POC(result_STPCA_MP_UCIDSA,DiscFea);
mean_POC_STPCA_MP_UCIDSA_disp = sprintf('%.2f%%', mean_POC_STPCA_MP_UCIDSA*100);
disp(['Mean POC of STPCA-MP on UCIDSA: ',mean_POC_STPCA_MP_UCIDSA_disp])

%%%%%%%%%%%% UCIHAR %%%%%%%%%%%%%%%%%
data_name = 'UCIHAR';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% obtain features with the greatest Between-class variance (BCV) of data
NumFS = 100;
[~,DiscFea] = BCV(data, label, tensor_type, NumFS);
% conduct STPCA-MP
para_special_STPCA_MP_UCIHAR = para_special_STPCA_MP;
para_special_STPCA_MP_UCIHAR.orientation = 2;
struct_STPCA_MP_UCIHAR = struct_STPCA_MP;
struct_STPCA_MP_UCIHAR.NumFS = NumFS;
struct_STPCA_MP_UCIHAR.para_special = para_special_STPCA_MP_UCIHAR;
[result_STPCA_MP_UCIHAR] = AlgExecution(struct_STPCA_MP_UCIHAR);
[mean_POC_STPCA_MP_UCIHAR,result_STPCA_MP_UCIHAR] = POC(result_STPCA_MP_UCIHAR,DiscFea);
mean_POC_STPCA_MP_UCIHAR_disp = sprintf('%.2f%%', mean_POC_STPCA_MP_UCIHAR*100);
disp(['Mean POC of STPCA-MP on UCIHAR: ',mean_POC_STPCA_MP_UCIHAR_disp])

%%%%%%%%%%%% UMCS %%%%%%%%%%%%%%%%%
data_name = 'UMCS';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% obtain features with the greatest Between-class variance (BCV) of data
NumFS = 15;
[~,DiscFea] = BCV(data, label, tensor_type, NumFS);
% conduct STPCA-MP
para_special_STPCA_MP_UMCS = para_special_STPCA_MP;
para_special_STPCA_MP_UMCS.orientation = 1;
struct_STPCA_MP_UMCS = struct_STPCA_MP;
struct_STPCA_MP_UMCS.NumFS = NumFS;
struct_STPCA_MP_UMCS.para_special = para_special_STPCA_MP_UMCS;
[result_STPCA_MP_UMCS] = AlgExecution(struct_STPCA_MP_UMCS);
[mean_POC_STPCA_MP_UMCS,result_STPCA_MP_UMCS] = POC(result_STPCA_MP_UMCS,DiscFea);
mean_POC_STPCA_MP_UMCS_disp = sprintf('%.2f%%', mean_POC_STPCA_MP_UMCS*100);
disp(['Mean POC of STPCA-MP on UMCS: ',mean_POC_STPCA_MP_UMCS_disp])

%%%%%%%%%%%%% ASCH %%%%%%%%%%%%%%%%%
data_name = 'ASCH';
data_path = fullfile('..', 'data', data_name);
config_Table7;
% obtain features with the greatest Between-class variance (BCV) of data
NumFS = 6;
[~,DiscFea] = BCV(data, label, tensor_type, NumFS);
% conduct STPCA-MP
para_special_STPCA_MP_ASCH = para_special_STPCA_MP;
para_special_STPCA_MP_ASCH.orientation = 2;
struct_STPCA_MP_ASCH = struct_STPCA_MP;
struct_STPCA_MP_ASCH.NumFS = NumFS;
struct_STPCA_MP_ASCH.para_special = para_special_STPCA_MP_ASCH;
[result_STPCA_MP_ASCH] = AlgExecution(struct_STPCA_MP_ASCH);
[mean_POC_STPCA_MP_ASCH,result_STPCA_MP_ASCH] = POC(result_STPCA_MP_ASCH,DiscFea);
mean_POC_STPCA_MP_ASCH_disp = sprintf('%.2f%%', mean_POC_STPCA_MP_ASCH*100);
disp(['Mean POC of STPCA-MP on ASCH: ',mean_POC_STPCA_MP_ASCH_disp])


%% create Table
create_table