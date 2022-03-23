clc;
close all;
clear all;

import org.opensim.modeling.*
%% Settings
cams_full_folder = 'data set';
% mkdir('..\opensim\experimental_data\');
output_folder = 'output data';

T = readtable('GRF_shift_times.csv');
T.Properties.RowNames = T.trial;
T = removevars(T, {'trial'});

n = 7;

temp=dir(fullfile(cams_full_folder,'*.csv'));

for i = 1:size(temp,1)

    csv_files{i} = [cams_full_folder '\' temp(i).name];
    
    csv_data = read_cams_full_csv(csv_files{i},-1,-1);
    temp1 = strsplit(csv_files{i},'\');
    filename = regexprep(temp(i,end).name,'_export_proc_data.csv','');
    
    % Finding the trial in GRF_shift_times.csv
    temp2 = strsplit(filename, '_');
    trial = [temp2(1) '_' temp2(end-1) '_' temp2(end)];
    trial = num2str(cell2mat(trial));
    
    temp_leg = strsplit(trial, '_');
    leg = temp_leg{1}(3);
    if leg == 'L'
        z_ang = 90;
        x_ang = -90;
        fp_seg = {'calcn_l','calcn_r','calcn_l','calcn_r','calcn_l','calcn_r','calcn_l','calcn_r','calcn_r'};
        shift_leg = 'R';
    else
        z_ang = -90;
        x_ang = -90;
        fp_seg = {'calcn_r','calcn_l','calcn_r','calcn_l','calcn_r','calcn_l','calcn_r','calcn_r','calcn_l'};
        shift_leg = 'L';
    end

    x_rot = [1,0,0; 0, cosd(x_ang), -sind(x_ang); 0, sind(x_ang), cosd(x_ang)]';
    z_rot = [cosd(z_ang),-sind(z_ang), 0; sind(z_ang), cosd(z_ang), 0; 0,0,1]';

    rot = z_rot*x_rot;
        
    HS_1 = T(trial, :).t_HS_1;
    HS_2 = T(trial, :).t_HS_2;
    FO_1 = T(trial, :).t_FO_1;
    FO_2 = T(trial, :).t_FO_2;
    
    csv_data.force_plates(9) = csv_data.force_plates(n);
    
    % Shifting COP
    COP_f = COP_shift(csv_data, FO_2, FO_1, HS_2, HS_1, n, shift_leg);
    csv_data.force_plates(9).COP_ifb(:,1:3) = COP_f;
    
    % Shifting F 
    F_f = F_shift(csv_data, FO_2, FO_1, HS_2, HS_1, n);
    csv_data.force_plates(9).F(:,1:3) = F_f;

    %% Write Moderated GRF 
    % Note: the order of gait events should be adapted based on the trial
    % being processed!
    fp_ind = [1,2,3,4,5,6,7,8,9];
    
    num_force_plates = length(fp_ind);
    grf_time =  csv_data.time;
    nGRFTime = length(grf_time);
    nPTs = length(1:nGRFTime);
    grf_time_vec = StdVectorDouble(nPTs);
    
    cc=0;
    for j = 1:nGRFTime
        grf_time_vec.set(cc,grf_time(j,1));
        cc=cc+1;
    end
    
    grf_mat = Matrix(nPTs,num_force_plates*9);
    
    cc=0;
    for t = 1:nGRFTime
        nFRCData = 0;

        for j = 1:length(fp_ind)
            
            force = csv_data.force_plates(fp_ind(j)).F(t,:)*rot;
            point = csv_data.force_plates(fp_ind(j)).COP_ifb(t,:)/1000*rot;
            torque = [0 0 csv_data.force_plates(fp_ind(j)).Tz_ifb(t)/1000]*rot;
            
            force(isnan(force))=0;
            point(isnan(point))=0;
            torque(isnan(torque))=0;

            for k = 1:3
                if (force(:,2) > 100)
                    grf_mat.set(cc,nFRCData,force(k));
                else
                    grf_mat.set(cc,nFRCData,0);
                end
                    nFRCData=nFRCData+1;
            end
            
            for k = 1:3
                if(force(:,2) > 100)
                    grf_mat.set(cc,nFRCData,point(k));
                else
                    grf_mat.set(cc,nFRCData,0);
                end
                nFRCData=nFRCData+1;
            end
            for k = 1:3
                if(force(:,2) > 100)
                    grf_mat.set(cc,nFRCData,torque(k));
                else
                    grf_mat.set(cc,nFRCData,0);
                end
                nFRCData=nFRCData+1;
            end            
        end
        cc=cc+1;
    end
   
    grf_names_vec = StdVectorString(num_force_plates*9);
    nFRCNames = 0;
    for j = 1:num_force_plates
        
        grf_names_vec.set((j-1)*9,[ 'ground_force_' int2str(fp_ind(j)) '_vx']);
        grf_names_vec.set((j-1)*9+1,['ground_force_' int2str(fp_ind(j)) '_vy']);
        grf_names_vec.set((j-1)*9+2,['ground_force_' int2str(fp_ind(j)) '_vz']);
        
        grf_names_vec.set((j-1)*9+3,['ground_force_' int2str(fp_ind(j)) '_px']);
        grf_names_vec.set((j-1)*9+4,['ground_force_' int2str(fp_ind(j)) '_py']);
        grf_names_vec.set((j-1)*9+5,['ground_force_' int2str(fp_ind(j)) '_pz']);
        
        grf_names_vec.set((j-1)*9+6,['ground_torque_' int2str(fp_ind(j)) '_x']);
        grf_names_vec.set((j-1)*9+7,['ground_torque_' int2str(fp_ind(j)) '_y']);
        grf_names_vec.set((j-1)*9+8,['ground_torque_' int2str(fp_ind(j)) '_z']);
        
        
    end
        
    grf_table = TimeSeriesTable(grf_time_vec, grf_mat, grf_names_vec);
    grf_table.addTableMetaDataString('nRows',int2str(nPTs));
    grf_table.addTableMetaDataString('nColumns',int2str(num_force_plates*9+1));
    grf_table.addTableMetaDataString('inDegrees','no');
    
    grf_mot_file = [filename '_grf.mot'];
    STOFileAdapter.write(grf_table, [output_folder grf_mot_file]);
    
 %% Write .xml ExternalLoads Files
    ext_loads = ExternalLoads();
    ext_loads.setDataFileName(grf_mot_file);    

    for j = 1:num_force_plates
        ext_frc = ExternalForce();
        ext_frc.setName(['FP' num2str(fp_ind(j))]);
        ext_frc.set_applied_to_body(fp_seg{j});
        ext_frc.set_force_expressed_in_body('ground');
        ext_frc.set_point_expressed_in_body('ground');
        ext_frc.set_force_identifier(['ground_force_' num2str(fp_ind(j)) '_v']);
        ext_frc.set_point_identifier(['ground_force_' num2str(fp_ind(j)) '_p']);
        ext_frc.set_torque_identifier(['ground_torque_' num2str(fp_ind(j)) '_']);
        ext_loads.cloneAndAppend(ext_frc);    
    end
    
    ext_load_file = [filename '_external_loads.xml'];
    ext_loads.print([output_folder ext_load_file]);
%     next_file = input('Enter...\n');
end