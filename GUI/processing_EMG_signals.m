%%% Function to process the data from the collaboraive robot to construct
%%% the pool of healthy synergies

function [W, H, nW, emg_trial, VAF, mean_emg, H_all, W_all, VAF2] = processing_KUKA_GUI_final(data, error, fs_emg, fs_acc)
%% Divide in EMG and ACC data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% EMG DATA %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the list of all column names from the dataTable
columnNames = data.Properties.VariableNames;

% Find columns starting with the names of the muscles
startsWithUpperTrapezius = startsWith(columnNames, 'UpperTrapezius');
startsWithPectoralisMajor = startsWith(columnNames, 'PectoralisMajor');
startsWithDeltoidAnterior = startsWith(columnNames, 'DeltoidAnterior');
startsWithDeltoidMedial = startsWith(columnNames, 'DeltoidMedial');
startsWithDeltoidPosterior = startsWith(columnNames, 'DeltoidPosterior');
startsWithBiceps = startsWith(columnNames, 'Biceps');
startsWithTriceps = startsWith(columnNames, 'Triceps');

% Combine the logical indices
selectedColumns = startsWithUpperTrapezius|startsWithPectoralisMajor|startsWithDeltoidAnterior|startsWithDeltoidMedial|startsWithDeltoidPosterior|startsWithBiceps|startsWithTriceps;
emg_data = data(:, selectedColumns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ACC DATA %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marker1 = startsWith(columnNames, 'Marker1');
marker2 = startsWith(columnNames, 'Marker2');
m1 = find(marker1); m2 = find(marker2);
idx = [m1:m1+3,m2:m2+3];
% acc_varnames = ["Var2","Var3","Var4","Var6","Var7","Var8"];
% acc_varnames = ["Var2","Var3","Var4","Var8","Var9","Var10"];
acc_data_pre = data(:,idx);

% Find rows with NaN values in any variable
rowsWithNaN = any(isnan(acc_data_pre{:,:}), 2);
% Remove rows with NaN values
acc_data = acc_data_pre(~rowsWithNaN, :);
t_acc = linspace(0,size(acc_data,1)/fs_acc,size(acc_data,1));
%% Preprocess (filtering & smoothing)

%%% Muscle activity envelope
t_emg = linspace(0,size(emg_data,1)/fs_emg,size(emg_data,1));
fprintf('max time emg %f \n',max(t_emg))
fprintf('max time acc %f \n\n',max(t_acc))

% ---- Upper trapezius envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'UpperT')));
[trapezius_env, ~, freq_tp, power_tp] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Pectoralis major envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'Pectoralis')));
[pectoralis_env, ~, freq_pecm, power_pecm] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Deltoid anterior envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'DeltoidAnt')));
[deltoid_ant_env, ~, freq_da, power_da] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Deltoid medial envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'DeltoidMed')));
[deltoid_med_env, ~, freq_dm, power_dm] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Deltoid posterior envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'DeltoidPost')));
[deltoid_pos_env, ~, freq_dp, power_dp] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Biceps envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'Biceps')));
[biceps_env, ~, freq_bic, power_bic] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
% ---- Triceps envelope
var_x = table2array(emg_data(:,contains(emg_data.Properties.VariableNames,'Triceps')));
[triceps_env, ~, freq_tri, power_tri] = methods_envelope_ul(var_x,fs_emg,'literature','spikes');
%% Normalize the envelopes in magnitude
% The moving mean value is calculated with a window of 96 points, which
% translates to ~50 ms (Ubeda et al., 2018).
val_norm = round(fs_emg * 0.05);
max_vals = [max(movmean(trapezius_env,val_norm)),max(movmean(pectoralis_env,val_norm)),...
    max(movmean(deltoid_ant_env,val_norm)),max(movmean(deltoid_med_env,val_norm)),max(movmean(deltoid_pos_env,val_norm)),...
    max(movmean(biceps_env,val_norm)),max(movmean(triceps_env,val_norm))];
power_vals = round([power_tp,power_pecm,power_da,power_dm,power_dp,power_bic,power_tri],4);

% ---- Upper trapezius envelope
if or(freq_tp <= 50, power_vals(1) < 0.0001)
    trapezius_norm = trapezius_env / 10;
elseif power_vals(1) == 0.0001
    trapezius_norm = trapezius_env / mean(max_vals);
else
    trapezius_norm = trapezius_env / max(movmean(trapezius_env,107));
end
% ---- Pectoralis major envelope
if or(freq_pecm <= 50, power_vals(2) < 0.0001)
    pectoralis_norm = pectoralis_env / 10;
elseif power_vals(2) == 0.0001    
    pectoralis_norm = pectoralis_env / mean(max_vals);
else
    pectoralis_norm = pectoralis_env / max(movmean(pectoralis_env,107));
end
% ---- Deltoid anterior envelope
if or(freq_da <= 50, power_vals(3) < 0.0001)
    deltoid_ant_norm = deltoid_ant_env / 10;
elseif power_vals(3) == 0.0001
    deltoid_ant_norm = deltoid_ant_env / mean(max_vals);
else
    deltoid_ant_norm = deltoid_ant_env / max(movmean(deltoid_ant_env,107));
end
% ---- Deltoid medial envelope
if or(freq_dm <= 50, power_vals(4) < 0.0001)
    deltoid_med_norm = deltoid_med_env / 10;
elseif power_vals(4) == 0.0001
    deltoid_med_norm = deltoid_med_env / mean(max_vals);
else
    deltoid_med_norm = deltoid_med_env / max(movmean(deltoid_med_env,107));
end
% ---- Deltoid posterior envelope
if or(freq_dp <= 50, power_vals(5) < 0.0001)
    deltoid_pos_norm = deltoid_pos_env / 10;
elseif power_vals(5) == 0.0001
    deltoid_pos_norm = deltoid_pos_env / mean(max_vals);
else
    deltoid_pos_norm = deltoid_pos_env / max(movmean(deltoid_pos_env,107));
end
% ---- Biceps envelope
if or(freq_bic <= 50, power_vals(6) < 0.0001)
    biceps_norm = biceps_env / 10;
elseif power_vals(6) == 0.0001
    biceps_norm = biceps_env / mean(max_vals);
else
    biceps_norm = biceps_env / max(movmean(biceps_env,107));
end
% ---- Triceps envelope
if or(freq_tri <= 50, power_vals(7) < 0.0001)
    triceps_norm = triceps_env / 10;
elseif power_vals(7) == 0.0001
    triceps_norm = triceps_env / mean(max_vals);
else
    triceps_norm = triceps_env / max(movmean(triceps_env,107));
end
%% Downsample the signals (if necessary for any processing)
factor = 3;
trapezius = downsample(trapezius_norm,factor);
pectoralis = downsample(pectoralis_norm,factor);
deltoid_ant = downsample(deltoid_ant_norm,factor);
deltoid_med = downsample(deltoid_med_norm,factor);
deltoid_pos = downsample(deltoid_pos_norm,factor);
biceps = downsample(biceps_norm,factor);
triceps = downsample(triceps_norm,factor);
%----- New time vector
t_emg_new = linspace(0,max(t_emg),length(biceps));
%% smooth acceleration
if fs_acc == 370
    var3 = table2array(acc_data(:,2)) + (0 - mean(table2array(acc_data(1:100,2))));
    var7 = table2array(acc_data(:,6)) + (0 - mean(table2array(acc_data(1:100,6))));
elseif or(fs_acc == 74,fs_acc == 148)
    var3 = table2array(acc_data(:,3)) + (0 - mean(table2array(acc_data(1:100,3))));
    var7 = table2array(acc_data(:,7)) + (0 - mean(table2array(acc_data(1:100,7))));
end
svar3 = smoothdata(var3);
svar7 = smoothdata(var7);
%% Identify start and stop labels
emg = [trapezius';pectoralis';deltoid_ant';deltoid_med';deltoid_pos';biceps';triceps'];
time = t_acc;

[pks,locs] = findpeaks(svar3,time,'MinPeakHeight',mean(svar3(svar3>0.005)),'MinPeakDistance',6,'MinPeakProminence',mean(svar3(svar3>0.03)));
[pks2,locs2] = findpeaks(svar7,time,'MinPeakHeight',mean(svar7(svar7>0.005)),'MinPeakDistance',6,'MinPeakProminence',mean(svar7(svar7>0.03)));

if error == 1
    locs(1) = [];   
    locs2(1) = [];
    pks(1) = [];
    pks2(1) = [];
elseif error == 2
    locs(9) = [];
    pks(9) = [];
elseif error ==3
    locs(end-1:end) = [];  
    pks(end-1:end) = [];
elseif error == 4
    locs(8) = [];  
    pks(8) = [];
    locs(6) = [];  
    pks(6) = [];
elseif error == 5
    locs(8) = [];
    pks(8) = [];
elseif error == 6
    locs(1) = [];
    pks(1) = [];
elseif error == 7
    locs(2) = [];
    pks(2) = [];
elseif error == 11
    pks = pks(4:end);
    locs = locs(4:end);
    pks2 = pks2(4:end);
    locs2 = locs2(4:end);
elseif error == 12
    locs(3) = [];
    pks(3) = []; 

end

% Delete the first 3 trials (assume as familiarization)

%----- Figure to visualize if the peaks are properly selected
% figure
% plot(time,svar3)
% hold on
% plot(time,svar7)
% hold on
% plot(locs,pks,'*')
% hold on
% plot(locs2,pks2,'*')
% legend('start','stop')
locs = locs(end-3:end);
pks = pks(end-3:end);
locs2 = locs2(end-3:end);
pks2 = pks2(end-3:end);
%% SYNERGIES
% 1) Calculate synergy per repetition and average
% 2) Average the reps and calculate synergy
%% 1) Calculate synergy per repetition and average
n_fs = 100;
fix_val = 5 * n_fs; % Approximately the task takes 5segs. Then 5 * 100 = 500Hz
t_fix = linspace(0,100,fix_val);
n = 3;
for i = 1:length(locs)
    % Find the index of the minimum difference
    [~, start_point] = min(abs(t_emg_new - locs(i)));
    [~, end_point] = min(abs(t_emg_new - locs2(i))); 
    t1 = linspace(0,100,length(t_emg_new(start_point:end_point)));
    xvals = emg(:,start_point:end_point);
    fix_val = (length(xvals)/fs_emg)*n_fs; % Number of points with 100Hz
    t_down = linspace(0,100,fix_val); % t vector for a new vector with fs=100hz
    emg_down = interp1(t1,xvals',t_down); % Downsample to 10ms - 100Hz
    x_emg = interp1(t_down,emg_down,t_fix);
    
    % ----- Calculate the NNMF
    % ----- Non-negative matrix
    rng(1)
    opt = statset('Display','final');
    [H2,W2] = nnmf(x_emg,n,'Replicates',100,'Options',opt);

    W(:,:,i) = W2';
    H(:,:,i) = H2';
    nW(:,:,i) = W2';  
    
    %----- Calculate the variance of muscle patterns
    % If one wants to identify VAF for a different n then change the input
    % in the nnmf function and run it again, is an iterative process.
    SST = sum(sum((x_emg' - mean(x_emg',2)).^2,2));
    SSE = sum(sum((x_emg' - W(:,:,i)*H(:,:,i)).^2,2));
    VAF(i) = 1 - (SSE/SST);
    
    % Matriz created to save the emg data divided per trial
    emg_trial(:,:,i) = x_emg;
end
mean_emg = mean(emg_trial,3);
rng(1)
[H_all, W_all] = nnmf(mean_emg,n,'Replicates',100,'Options',opt);
SST = sum(sum((mean_emg - mean(mean_emg)).^2));
SSE = sum(sum((mean_emg' - W_all'*H_all').^2));
VAF2 = 1 - (SSE/SST);
% Dimensions W: 7 muscles X n synergies X no reps(4)
% Dimensions H: n ac X 1000 (fix len) X no reps(4)
end