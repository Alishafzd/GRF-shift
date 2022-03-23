function data = read_cams_full_csv(file,start_time,stop_time)
%% read_cams_sample_csv
%==========================================================================
% Author: Colin Smith - Hamed Hosseini
%
% Read the .csv files provided in the full CAMS-Knee datasets and return
% a MATLAB structure containing the data.
%
% Set start_time or stop_time to use the first and last times listed in the
% .csv file.
%==========================================================================

%% Open File
fid=fopen(file);
if fid==-1
    printf('File not found: %s',file)
end

nEvent = 0;

%% Read Header
for i = 1:25
    line = fgetl(fid);
    
    if(i==23)
        data_type = strsplit(line,',','CollapseDelimiters',false);
        
        marker_start = find(strcmp(data_type,'TRAJECTORIES'));
        fp(1).start = find(strcmp(data_type,'FORCE PLATE 1'));
        fp(2).start = find(strcmp(data_type,'FORCE PLATE 2'));
        fp(3).start = find(strcmp(data_type,'FORCE PLATE 3'));
        fp(4).start = find(strcmp(data_type,'FORCE PLATE 4'));
        fp(5).start = find(strcmp(data_type,'FORCE PLATE 5'));
        fp(6).start = find(strcmp(data_type,'FORCE PLATE 6'));
        fp(7).start = find(strcmp(data_type,'FORCE PLATE 7'));
        fp(8).start = find(strcmp(data_type,'FORCE PLATE 8'));
        emg_start = find(strcmp(data_type,'EMG DATA'));
        knee_load_start = find(strcmp(data_type,'In-vivo knee prosthesis data'));
        tibia_trans_start = find(strcmp(data_type,'Transformation matrix tibia to lab CS'));
        femur_trans_start = find(strcmp(data_type,'Transformation matrix femur to lab CS'));
        
        marker_end = fp(1).start-1;
        fp(1).end = fp(2).start-1;
        fp(2).end = fp(3).start-1;
        fp(3).end = fp(4).start-1;
        fp(4).end = fp(5).start-1;
        fp(5).end = fp(6).start-1;
        fp(6).end = fp(7).start-1;
        fp(7).end = fp(8).start-1;
        fp(8).end = emg_start-1;
        emg_end = knee_load_start-1;
        knee_load_end = tibia_trans_start-1;
        tibia_trans_end = femur_trans_start-1;
        femur_trans_end = length(data_type);
        
    elseif(i==24)
        labels = strsplit(line,',','CollapseDelimiters',false);
    elseif(i==25)
        units = strsplit(line,',','CollapseDelimiters',false);
    end
    
    %Read Gait Events
    if (i>9 && i < 23)
        if isempty(line); continue; end;
        nEvent = nEvent+1;
        split_line = strsplit(line,',');
        events(nEvent).name = split_line{1};
        events(nEvent).time = split_line{2};
    end
        
end

raw_data = csvread(file,26,0);
nSamples = size(raw_data,1);



%% Set start and stop times
if start_time == -1 
    start_frame = 1;
else
    start_frame = interp1(raw_data(:,1),1:length(raw_data(:,1)),start_time,'nearest');
end

if stop_time == -1 
    stop_frame = length(raw_data(:,1));
else
    stop_frame = interp1(raw_data(:,1),1:length(raw_data(:,1)),stop_time,'nearest');
end

data.time = raw_data(start_frame:stop_frame,1);

%% Marker Trajectories
nMarkers=0;

marker_sample = zeros(nSamples,1);

%Find marker samples in time 
for i = start_frame:stop_frame
    if (raw_data(i,3) == floor(raw_data(i,3)))
        marker_sample(i) = 1;
    end
end

nMarkerSample = nnz(marker_sample(start_frame:stop_frame));

for i = marker_start:3:marker_end
    if ~isempty(labels{i})
        nMarkers = nMarkers+1;
        marker_name = regexprep(labels{i},{'_X','_Y','_Z'},'');
        data.marker_trajectories.marker_names{nMarkers} = marker_name;
        data.marker_trajectories.(marker_name) = nan(nMarkerSample,3);
    end
    
    
    marker_comp = units{i}(1);
    
    nMarkerSample = 1;
    for j = 1: nSamples
        
        if (j < start_frame)
            continue;
        elseif(j > stop_frame)
            break;
        end
            
        if (marker_sample(j) == 0)
            continue;
        end
        
        
        
        if nMarkers == 1 
            data.marker_trajectories.time(nMarkerSample) = raw_data(j,1);
        end        

            data.marker_trajectories.(marker_name)(nMarkerSample,1) = raw_data(j,i);
            data.marker_trajectories.(marker_name)(nMarkerSample,2) = raw_data(j,i+1);
            data.marker_trajectories.(marker_name)(nMarkerSample,3) = raw_data(j,i+2);
            
        nMarkerSample=nMarkerSample+1;
    end
    
end
data.marker_trajectories.nMarkers = nMarkers;

%% Force Plates
for k = 1:8
    
    data.force_plates(k).ref(:,1) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Ref_X_lab_' int2str(k)])));
    data.force_plates(k).ref(:,2) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Ref_Y_lab_' int2str(k)])));
    data.force_plates(k).ref(:,3) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Ref_Z_lab_' int2str(k)])));
    
    data.force_plates(k).F(:,1) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Fx_lab_' int2str(k)])));
    data.force_plates(k).F(:,2) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Fy_lab_' int2str(k)])));
    data.force_plates(k).F(:,3) = raw_data(start_frame:stop_frame,find(strcmp(labels,['Fz_lab_' int2str(k)])));
    
    data.force_plates(k).COP(:,1) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPx_lab_' int2str(k)])));
    data.force_plates(k).COP(:,2) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPy_lab_' int2str(k)])));
    data.force_plates(k).COP(:,3) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPz_lab_' int2str(k)])));
    
    data.force_plates(k).Tz = raw_data(start_frame:stop_frame,find(strcmp(labels,['Tz_lab_' int2str(k)])));
    
    data.force_plates(k).COP_ifb(:,1) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPx_lab_IfB_' int2str(k) ])));
    data.force_plates(k).COP_ifb(:,2) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPy_lab_IfB_' int2str(k) ])));    
    data.force_plates(k).COP_ifb(:,3) = raw_data(start_frame:stop_frame,find(strcmp(labels,['COPz_lab_IfB_' int2str(k) ])));
    
    data.force_plates(k).Tz_ifb = raw_data(start_frame:stop_frame,find(strcmp(labels,['Tz_lab_IfB_' int2str(k)])));
    
end


%% EMG
    data.emg.left.rectus_femoris = raw_data(start_frame:stop_frame,strcmp(labels,'RectFemL_01'));
    data.emg.left.vastus_medialis = raw_data(start_frame:stop_frame,strcmp(labels,'VastusMedL_02'));
    data.emg.left.vastus_lateralis = raw_data(start_frame:stop_frame,strcmp(labels,'VastusLatL_03'));
    data.emg.left.tibialis_anterior = raw_data(start_frame:stop_frame,strcmp(labels,'TibAntL_04'));
    data.emg.left.medial_hamstrings = raw_data(start_frame:stop_frame,strcmp(labels,'HamMedL_05'));
    data.emg.left.lateral_hamstrings = raw_data(start_frame:stop_frame,strcmp(labels,'HamLatL_06'));
    data.emg.left.medial_gastrocnemius = raw_data(start_frame:stop_frame,strcmp(labels,'GastroMedL_07'));
    data.emg.left.lateral_gastrocnemius = raw_data(start_frame:stop_frame,strcmp(labels,'GastroLatL_08'));
    data.emg.right.rectus_femoris = raw_data(start_frame:stop_frame,strcmp(labels,'RectFemR_09'));
    data.emg.right.vastus_medialis = raw_data(start_frame:stop_frame,strcmp(labels,'VastusMedR_10'));
    data.emg.right.vastus_lateralis = raw_data(start_frame:stop_frame,strcmp(labels,'VastusLatR_11'));
    data.emg.right.tibialis_anterior = raw_data(start_frame:stop_frame,strcmp(labels,'TibAntR_12'));
    data.emg.right.medial_hamstrings = raw_data(start_frame:stop_frame,strcmp(labels,'HamMedR_13'));
    data.emg.right.lateral_hamstrings = raw_data(start_frame:stop_frame,strcmp(labels,'HamLatR_14'));
    data.emg.right.medial_gastrocnemius = raw_data(start_frame:stop_frame,strcmp(labels,'GastroMedR_15'));
    data.emg.right.lateral_gastrocnemius = raw_data(start_frame:stop_frame,strcmp(labels,'GastroLatR_16'));

%% Knee Loads
knee_load_sample = zeros(nSamples,1);

%Find knee load samples in time 
for t = start_frame:stop_frame
    if (~isnan(raw_data(t,knee_load_start)))
        knee_load_sample(t) = 1;
    end
end

nKneeLoadSample = nnz(knee_load_sample(start_frame:stop_frame));
data.knee_loads.time = zeros(nKneeLoadSample,1);

for i = knee_load_start:knee_load_end

    knee_load_comp = labels{i};
    data.knee_loads.(knee_load_comp) = zeros(nKneeLoadSample,1);
    
    kneeSampleInd=1;
    for j = 1: nSamples
        if (j < start_frame)
            continue;
        elseif(j > stop_frame)
            break;
        end

        if (knee_load_sample(j) == 0)
            continue;
        end

        if i == knee_load_start
            data.knee_loads.time(kneeSampleInd) = raw_data(j,1);
        end

        data.knee_loads.(knee_load_comp)(kneeSampleInd) = raw_data(j,i);
        kneeSampleInd=kneeSampleInd+1;
    end    
end

%% Fluoroscopy Knee Kinematics
fluoro_sample_col = find(strcmp(labels,'SYNC_fluoro'));
fluoro_sample = raw_data(:,fluoro_sample_col);
nFluoroSamples = nnz(fluoro_sample(start_frame:stop_frame));
data.knee_kinematics.time = zeros(nFluoroSamples,1);

for i = tibia_trans_start:tibia_trans_end
    trans_label = strrep(labels{i},'-','_');
    
    data.knee_kinematics.(trans_label) = zeros(nFluoroSamples,1);
    fluoroSampleInd=1;
    
    for j = 1: nSamples
        if (j < start_frame)
            continue;
        elseif(j > stop_frame)
            break;
        end

        if (fluoro_sample(j) == 0)
            continue;
        end
        
        if i == tibia_trans_start
            data.knee_kinematics.time(fluoroSampleInd) = raw_data(j,1);
        end
        
        data.knee_kinematics.(trans_label)(fluoroSampleInd) = raw_data(j,i);
        fluoroSampleInd=fluoroSampleInd+1;
    end
end

for i = femur_trans_start:femur_trans_end
    trans_label = strrep(labels{i},'-','_');
    
    data.knee_kinematics.(trans_label) = zeros(nFluoroSamples,1);
    fluoroSampleInd=1;
    
    for j = 1: nSamples
        if (j < start_frame)
            continue;
        elseif(j > stop_frame)
            break;
        end

        if (fluoro_sample(j) == 0)
            continue;
        end
        
        data.knee_kinematics.(trans_label)(fluoroSampleInd) = raw_data(j,i);
        fluoroSampleInd=fluoroSampleInd+1;
    end
end
fclose(fid);