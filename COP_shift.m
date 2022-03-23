function [COP_f] = COP_shift(csv_data,t_i_end,t_f_end,t_i_start,t_f_start,n,leg)

toe = fillgaps(csv_data.marker_trajectories.([leg 'TTO']));
heel = fillgaps(csv_data.marker_trajectories.([leg 'THL']));
ankle = fillgaps(csv_data.marker_trajectories.([leg 'TVM']));
COP = fillgaps(csv_data.force_plates(n).COP_ifb);

% Resampling markers to COP length
toe_resample = zeros(length(COP), 3);
for j = 1:3
    toe_resample(:, j) = MinaSize(toe(:, j), length(COP));
end

heel_resample = zeros(length(COP), 3);
for j = 1:3
    heel_resample(:, j) = MinaSize(heel(:, j), length(COP));
end

ankle_resample = zeros(length(COP), 3);
for j = 1:3
    ankle_resample(:, j) = MinaSize(ankle(:, j), length(COP));
end

% Changing coordinate system to local 
localCOP = COP_local(heel_resample,toe_resample,ankle_resample,COP);

% Resampling and shifting the COP in local coordinate system
m = round(2000*(t_f_end-t_f_start))/round(2000*(t_i_end-t_i_start));
localCOP_resample = zeros(round(m*length(localCOP)), 3);
for j = 1:3
    localCOP_resample(:, j) = MinaSize(localCOP(:, j), round(m*length(localCOP)));
end

t_i_end_m = m*t_i_end;
shift = round(2000*(t_f_end-t_i_end_m));

t = 1:1:length(localCOP_resample);
t_resample = t + shift*ones(1,length(t));

localCOP2 = localCOP_resample(find(t_resample > 0):1:length(t_resample),1:3);
[b, a] = butter(4, 40/2000, 'low');
localCOP2 = filter(b, a, localCOP2);

% Changing COP local coordinate system to global
COP_temp = COP_global(heel_resample,toe_resample,ankle_resample,localCOP2);

% Chaning dimension to standard
h = find(csv_data.force_plates(n).COP_ifb(:,3) > 0);
delta_h = h(end)-h(1);
j = round(t_f_end*2000)-delta_h*m;

COP_f = zeros(length(csv_data.force_plates(n).COP_ifb),3);
COP_f(round(j):round(t_f_end*2000),:) = COP_temp(round(j):round(t_f_end*2000),:);

end
