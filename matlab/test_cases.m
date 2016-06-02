clear;


file_path = '/Users/zahedi/bwSyncAndShare/Article_QuantMorphComp/Model_Results_Data/';
suffix    = '_noPert_h7cm';

% files(1).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_motor'                    suffix];
% files(2).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_linFv_constFl_FFB' suffix];
% files(3).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_HillFv_HillFl_FFB' suffix];

files(1).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_motor'                     suffix];
files(2).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_linFv_constFl_FFB'  suffix];
files(3).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_linFv_linFl_FFB'    suffix];
files(4).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_HillFv_constFl_FFB' suffix];
files(5).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_HillFv_HillFl_FFB'  suffix];
files(6).name = [file_path '/NoPerturbationSameHoppingHeight/' 'results_muscle_MTC_FFB'            suffix];


% 1. Step: Extract min/max values for each column that will be used
p_min  = 0; p_max = 0;
ac_min = 0; ac_max = 0;
v_min  = 0; v_max = 0;
a_min  = 0; a_max = 0;
s_min  = 0; s_max = 0;

domain_string = '';


for file_index = 1:length(files)
    load(files(file_index).name,'SimData');
    [pathstr,name,ext] = fileparts(files(file_index).name);
    files(file_index).short_name = name;
    
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
        
        a_min = min(action);
        a_max = max(action);
        
        s_min = min(muscle_sensor);
        s_max = max(muscle_sensor);
        
        ac_min = min(accelaration);
        ac_max = max(accelaration);
    else
        p_min = min(p_min, min(position));
        p_max = max(p_max, max(position));
        
        v_min = min(v_min, min(velocity));
        v_max = max(v_max, max(velocity));
        
        a_min = min(a_min, min(action));
        a_max = max(a_max, max(action));
        
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
    sprintf('  Actuator signal: %f, %f\n', a_min, a_max)];

%%
for bins = [30]
    
    % 2. Step: Discretise data
    w_bins = bins;
    a_bins = bins;
    s_bins = bins;
    
    tic
    
    %fprintf('Binning the data with |W| = %d, |A| = %d, |S| = %d\n', w_bins, a_bins, s_bins);
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        % fprintf('Discretising file %s\n', name);
        load(files(file_index).name, 'SimData');
        
        position      = SimData(:,2);
        velocity      = SimData(:,3);
        accelaration  = SimData(:,4);
        action        = SimData(:,10);
        muscle_sensor = SimData(:,5);
        
        files(file_index).d_position     = discretiseMatrix(position, p_min, p_max, w_bins);
        files(file_index).d_velocity     = discretiseMatrix(velocity, v_min, v_max, w_bins);
        files(file_index).d_accelaration = discretiseMatrix(accelaration, ac_min, ac_max, w_bins);
        
        files(file_index).a              = discretiseMatrix(action, a_min, a_max, a_bins);
        files(file_index).w              = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity, files(file_index).d_accelaration]);
        
        if isempty(strfind(files(file_index).name,'muscle'))
            % fprintf('%s is not a muscle\n', name);
            files(file_index).s          = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity]);
        else
            % fprintf('%s is a muscle\n', name);
            files(file_index).s          = discretiseMatrix(muscle_sensor, s_min, s_max, s_bins);
            
        end
    end
    
    %% MC_W Test case
    
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        fprintf('Working on file %s\n', name);
        w2 = files(file_index).w(2:end,:);
        w1 = files(file_index).w(1:end-1,:);
        a1 = files(file_index).a(1:end-1,:);
        s1 = files(file_index).s(1:end-1,:);
        
        %files(file_index).mcmi = MC_MI_dynamic(w2w1a1s1);
        files(file_index).mcw   = MC_W(w2, w1, a1);
        files(file_index).mcw1  = MC_W1(w2, w1, a1);
        %files(file_index).mcc  = MC_C_dynamic(wsa);
        
        if abs(files(file_index).mcw - files(file_index).mcw1) > 0.000001
            fprintf('Error in MC_W for %d bins in file %s\n', bins, files(file_index).short_name);
        else
            fprintf('MC_W for %d bins in file %s OK error is %f\n', bins, files(file_index).short_name,abs(files(file_index).mcw - files(file_index).mcw1));
        end
        
        files(file_index).mcwd  = MC_W_dynamic(w2, w1, a1);
        files(file_index).mcw1d = MC_W1_dynamic(w2, w1, a1);
        
        s = mean((files(file_index).mcwd - files(file_index).mcw1d).^2);
        
        if s > 0.000001
            fprintf('Error in MC_W_dynamic error for %d bins in file %s = %f\n', bins, files(file_index).short_name, s);
        else
            fprintf('MC_W_dynamic on %d bins in file %s = %f OK\n', bins, files(file_index).short_name, s);
        end
    end
    
    
%     fprintf('MC_W\n');
%     for file_index = 1:length(files)
%         fprintf('File %s\n', files(file_index).short_name);
%         fprintf('  MC_W:  %f\n',   files(file_index).mcw);
%         fprintf('  MC_W1: %f\n',   files(file_index).mcw1);
%     end
    
    %% MC_MI Test case
    
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        w2 = files(file_index).w(2:end,:);
        w1 = files(file_index).w(1:end-1,:);
        a1 = files(file_index).a(1:end-1,:);
        s1 = files(file_index).s(1:end-1,:);
        
        %files(file_index).mcmi = MC_MI_dynamic(w2w1a1s1);
        files(file_index).mcmi   = MC_MI(w2, w1, s1, a1);
        files(file_index).mcmi1  = MC_MI1(w2, w1, s1, a1);
        files(file_index).mcmi2  = MC_MI2(w2, w1, s1, a1);
        %files(file_index).mcc  = MC_C_dynamic(wsa);
        
        if abs(files(file_index).mcmi - files(file_index).mcmi1) > 0.000001
            fprintf('Error in MC_MI for %d bins in file %s\n', bins, files(file_index).short_name);
        end
        if abs(files(file_index).mcmi - files(file_index).mcmi2) > 0.000001
            fprintf('Error in MC_MI for %d bins in file %s\n', bins, files(file_index).short_name);
        end
        
        files(file_index).mcmid = MC_MI_dynamic(w2, w1, s1, a1);
        files(file_index).mcmi1d = MC_MI1_dynamic(w2, w1, s1, a1);
        
        s = mean((files(file_index).mcmid - files(file_index).mcmi1d).^2);
        
        if s > 0.000001
            fprintf('Error in MC_MI_dynamic error for %d bins in file %s = %f\n', bins, files(file_index).short_name, s);
        else
            fprintf('MC_MI_dynamic in %d bins in file %s = %f OK\n', bins, files(file_index).short_name);
        end
    end
    
    
%     fprintf('MC_MI\n');
%     for file_index = 1:length(files)
%         fprintf('File %s\n', files(file_index).short_name);
%         fprintf('  MC_MI:  %f\n',   files(file_index).mcmi);
%         fprintf('  MC_MI1: %f\n',   files(file_index).mcmi1);
%         fprintf('  MC_MI2: %f\n',   files(file_index).mcmi2);
%     end
    
    
    %% MC_C Test case
    
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        w2 = files(file_index).w(2:end,:);
        w1 = files(file_index).w(1:end-1,:);
        a1 = files(file_index).a(1:end-1,:);
        s1 = files(file_index).s(1:end-1,:);
        
        files(file_index).mcc  = MC_C(w2, w1, s1, a1);
        files(file_index).mcc1 = MC_C1(w2, w1, s1, a1);
        
        if abs(files(file_index).mcc - files(file_index).mcc1) > 0.000001
            fprintf('Error in MC_C for %d bins in file %s\n', bins, files(file_index).short_name);
        else
            fprintf('MC_C for %d bins in file %s OK\n', bins, files(file_index).short_name);
        end
        
        files(file_index).mccd = MC_C_dynamic(w2, w1, s1, a1);
        files(file_index).mcc1d = MC_C1_dynamic(w2, w1, s1, a1);
        
        s = mean((files(file_index).mccd - files(file_index).mcc1d).^2);
        
        if s > 0.000001
            fprintf('Error in MC_C_dynamic error for %d bins in file %s = %f\n', bins, files(file_index).short_name, s);
        else
            fprintf('MC_C_dynamic on %d bins in file %s = %f OK\n', bins, files(file_index).short_name);
        end
    end
    
%     fprintf('MC_C\n');
%     for file_index = 1:length(files)
%         fprintf('File %s\n', files(file_index).short_name);
%         fprintf('  MC_C:  %f\n',   files(file_index).mcc);
%         fprintf('  MC_C1: %f\n',   files(file_index).mcc1);
%     end
end