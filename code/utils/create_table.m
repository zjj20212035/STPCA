%% summarize the result
best_acc_PIE = result_STPCA_MP_PIE.best_acc_disp;
best_nmi_PIE = result_STPCA_MP_PIE.best_nmi_disp;
best_ss_PIE = result_STPCA_MP_PIE.best_ss_disp;
time_PIE = result_STPCA_MP_PIE.time_disp;

best_acc_COIL20 = result_STPCA_MP_COIL20.best_acc_disp;
best_nmi_COIL20 = result_STPCA_MP_COIL20.best_nmi_disp;
best_ss_COIL20 = result_STPCA_MP_COIL20.best_ss_disp;
time_COIL20 = result_STPCA_MP_COIL20.time_disp;

best_acc_Imm40 = result_STPCA_MP_Imm40.best_acc_disp;
best_nmi_Imm40 = result_STPCA_MP_Imm40.best_nmi_disp;
best_ss_Imm40 = result_STPCA_MP_Imm40.best_ss_disp;
time_Imm40 = result_STPCA_MP_Imm40.time_disp;

best_acc_JAFFE = result_STPCA_MP_JAFFE.best_acc_disp;
best_nmi_JAFFE = result_STPCA_MP_JAFFE.best_nmi_disp;
best_ss_JAFFE = result_STPCA_MP_JAFFE.best_ss_disp;
time_JAFFE = result_STPCA_MP_JAFFE.time_disp;

best_acc_USPS = result_STPCA_MP_USPS.best_acc_disp;
best_nmi_USPS = result_STPCA_MP_USPS.best_nmi_disp;
best_ss_USPS = result_STPCA_MP_USPS.best_ss_disp;
time_USPS = result_STPCA_MP_USPS.time_disp;

best_acc_BreastMNIST = result_STPCA_MP_BreastMNIST.best_acc_disp;
best_nmi_BreastMNIST = result_STPCA_MP_BreastMNIST.best_nmi_disp;
best_ss_BreastMNIST = result_STPCA_MP_BreastMNIST.best_ss_disp;
time_BreastMNIST = result_STPCA_MP_BreastMNIST.time_disp;

time_UCIDSA = result_STPCA_MP_UCIDSA.time_disp;
time_UCIHAR = result_STPCA_MP_UCIHAR.time_disp;
time_UMCS = result_STPCA_MP_UMCS.time_disp;
time_ASCH = result_STPCA_MP_ASCH.time_disp;

%% print results
disp('The best result of STPCA-MP on tube-wise datasets: ');
headers_1 = {'Metric', 'PIE', 'COIL20', 'JAFFE', 'BreastMNIST', 'Imm40', 'USPS'};
clustering_metrics = {
    'ACC', best_acc_PIE, best_acc_COIL20, best_acc_JAFFE, best_acc_BreastMNIST, best_acc_Imm40, best_acc_USPS;
    'NMI', best_nmi_PIE, best_nmi_COIL20, best_nmi_JAFFE, best_nmi_BreastMNIST, best_nmi_Imm40, best_nmi_USPS; 
    'SS', best_ss_PIE, best_ss_COIL20, best_ss_JAFFE, best_ss_BreastMNIST, best_ss_Imm40, best_ss_USPS; 
    'TIME', time_PIE, time_COIL20, time_JAFFE, time_BreastMNIST, time_Imm40, time_USPS;            
};

fprintf('%-10s %-10s %-10s %-10s %-12s %-10s %-10s\n', headers_1{:});
fprintf('%-10s %-10s %-10s %-10s %-12s %-10s %-10s\n', '----------', '----------', '----------', '----------', '------------', '----------', '----------');

for i = 1:size(clustering_metrics, 1)
    metric_name = clustering_metrics{i, 1};
    values = clustering_metrics(i, 2:end);
    switch metric_name
        case {'ACC', 'NMI','SS'}
            formatted_values = values;
        case 'TIME'
            formatted_values = cellfun(@(x) sprintf('%.4fs', x), values, 'UniformOutput', false);
    end
    
    fprintf('%-10s ', metric_name);
    fprintf('%-10s %-10s %-10s %-12s %-10s %-10s\n', formatted_values{:});
end


disp('The best result of STPCA-MP on slice-wise datasets: ');
headers_2 = {'Metric', 'UCIDSA', 'UCIHAR', 'UMCS', 'ASCH'};
POC_metrics = {
    'POC', mean_POC_STPCA_MP_UCIDSA, mean_POC_STPCA_MP_UCIHAR, mean_POC_STPCA_MP_UMCS, mean_POC_STPCA_MP_ASCH;
    'TIME', time_UCIDSA, time_UCIHAR, time_UMCS, time_ASCH
};

fprintf('%-10s %-10s %-10s %-10s %-10s\n', headers_2{:});
fprintf('%-10s %-10s %-10s %-10s %-10s\n', ...
    '----------', '----------', '----------', '----------', '----------');

for i = 1:size(POC_metrics, 1)
    metric_name = POC_metrics{i, 1};
    values = POC_metrics(i, 2:end);
    
    switch metric_name
        case 'POC'
            formatted_values = cellfun(@(x) sprintf('%.2f%%', x*100), values, 'UniformOutput', false);
        case 'TIME'
            formatted_values = cellfun(@(x) sprintf('%.4fs', x), values, 'UniformOutput', false);
    end
   
    fprintf('%-10s %-10s %-10s %-10s %-10s\n', ...
        metric_name, formatted_values{:});
end