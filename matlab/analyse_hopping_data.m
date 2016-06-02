%% Initialization

clear;
bins = 300;

files(1).name = 'data/DCMot_ModelOutput';
files(2).name = 'data/MusLin_ModelOutput';
files(3).name = 'data/MusFib_ModelOutput';

%% 1. Step: Extract min/max values for each column that will be used
p_min  = 0; p_max = 0;
ac_min = 0; ac_max = 0;
v_min  = 0; v_max = 0;
a_min  = 0; a_max = 1.0;
s_min  = 0; s_max = 0;

domain_string = '';


for file_index = 1:length(files)
    load(files(file_index).name,'SimData');
    
    position      = SimData(:,2);
    velocity      = SimData(:,3);
    accelaration  = SimData(:,4);
    action        = SimData(:,10);
    muscle_sensor = SimData(:,5);
    
    files(file_index).position = position;
    
    if file_index == 1
        p_min = min(position);
        p_max = max(position);
        
        v_min = min(velocity);
        v_max = max(velocity);
        
        s_min = min(muscle_sensor);
        s_max = max(muscle_sensor);
        
        ac_min = min(accelaration);
        ac_max = max(accelaration);
    else
        p_min = min(p_min, min(position));
        p_max = max(p_max, max(position));
        
        v_min = min(v_min, min(velocity));
        v_max = max(v_max, max(velocity));
        
        ac_min = min(ac_min, min(accelaration));
        ac_max = max(ac_max, max(accelaration));
        
        s_min = min(s_min, min(muscle_sensor));
        s_max = max(s_max, max(muscle_sensor));
    end
end

domain_string = [domain_string 'Domains\n' ...
    sprintf('  Position:        %f, %f\n', p_min, p_max) ...
    sprintf('  Velocity:        %f, %f\n', v_min, v_max) ...
    sprintf('  Accelaration:    %f, %f\n', ac_min, ac_max) ...
    sprintf('  Muscle sensor:   %f, %f\n', s_min, s_max) ...
    sprintf('  Actuator signal: %f, %f\n', a_min, a_max)];



%%

fprintf('Bins = %d\n', bins);
w_bins = bins;
a_bins = bins;
s_bins = bins;

tic
for file_index = 1:length(files)
    [pathstr,name,ext] = fileparts(files(file_index).name);
    load(files(file_index).name, 'SimData');
    
    position      = SimData(:,2);
    velocity      = SimData(:,3);
    accelaration  = SimData(:,4);
    action        = SimData(:,10);
    muscle_sensor = SimData(:,5);
    
    if isempty(strfind(files(file_index).name,'muscle'))
        action = (action - min(action)) / (max(action) - min(action));
    end
    
    files(file_index).d_position     = discretiseMatrix(position, p_min, p_max, w_bins);
    files(file_index).d_velocity     = discretiseMatrix(velocity, v_min, v_max, w_bins);
    files(file_index).d_accelaration = discretiseMatrix(accelaration, ac_min, ac_max, w_bins);
    files(file_index).position       = position;
    files(file_index).velocity       = velocity;
    files(file_index).accelaration   = accelaration;
    files(file_index).action         = action;
    files(file_index).msensor        = muscle_sensor;
    
    files(file_index).a              = discretiseMatrix(action, a_min, a_max, a_bins);
    files(file_index).w              = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity, files(file_index).d_accelaration]);
    
    if isempty(strfind(files(file_index).name,'Mus'))
        fprintf('%s is not a muscle\n', name);
        files(file_index).s          = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity]);
    else
        fprintf('%s is a muscle\n', name);
        files(file_index).s          = discretiseMatrix(muscle_sensor, s_min, s_max, s_bins);
    end
    
end

%% calculate the measures

%fprintf('Starting calculations\n')
for file_index = 1:length(files)
    [pathstr,name,ext] = fileparts(files(file_index).name);
    files(file_index).short_name = name;
    
    fprintf('Working on file %s\n', name);
    w2 = files(file_index).w(2:end,:);
    w1 = files(file_index).w(1:end-1,:);
    a1 = files(file_index).a(1:end-1,:);
    s1 = files(file_index).s(1:end-1,:);
    
    fprintf('Calculating MC_W\n');
    tic
    files(file_index).mcw  = MC_W(w2, w1, a1);
    fprintf('result %f\n', files(file_index).mcw);
    mcwd = MC_W_dynamic(w2, w1, a1);
    files(file_index).mcd  = mcwd;
    fprintf('check %f\n', files(file_index).mcw - mean(mcwd));
    toc
    
    fprintf('Calculating MC_MI\n');
    tic
    files(file_index).mcmi = MC_MI(w2, w1, s1, a1);
    fprintf('result %f\n', files(file_index).mcmi);
    mcmid = MC_MI_dynamic(w2, w1, s1, a1);
    files(file_index).mcd  = [files(file_index).mcd mcmid];
    fprintf('check %f\n', files(file_index).mcmi - mean(mcmid));
    toc
    
end

%% write the data to csv files
csvwrite(sprintf('data/mc_data_dyn_dcmot_%d.csv', bins), files(1).mcd);
csvwrite(sprintf('data/mc_data_dyn_lin_%d.csv',   bins), files(2).mcd);
csvwrite(sprintf('data/mc_data_dyn_fib_%d.csv',   bins), files(3).mcd);

filename = sprintf('data/results_w%d_a%d_s%d.txt', w_bins, a_bins, s_bins);
fprintf('Writing results to %s\n', filename);

fileID = fopen(filename,'w');
fprintf(fileID, domain_string);
fprintf(fileID, 'Bins %d\n', bins);
for file_index = 1:length(files)
    [pathstr,name,ext] = fileparts(files(file_index).name);
    fprintf(fileID, 'Filename %s\n', name);
    fprintf(fileID, '  MC_W:  %f\n', files(file_index).mcw);
    fprintf(fileID, '  MC_MI: %f\n', files(file_index).mcmi);
end
fclose(fileID);
fprintf('done.\n')



%% Writing hopping data

for file_index = 1:length(files)
    [pathstr,name,ext] = fileparts(files(file_index).name);
    % fprintf('Discretising file %s\n', name);
    load(files(file_index).name, 'SimData');
    
    position      = SimData(:,2);
    velocity      = SimData(:,3);
    accelaration  = SimData(:,4);
    action        = SimData(:,10);
    
    r = [];
    r = [r position];
    r = [r velocity];
    r = [r accelaration];
    r = [r action];
    
    if isempty(strfind(files(file_index).name,'muscle')) == 0
        muscle_sensor = SimData(:,5);
        r = [r muscle_sensor];
    end
    
    fprintf('writing to data/hopping_data_%s.csv\n', name);
    csvwrite(sprintf('data/hopping_data_%s.csv', name), r);
end