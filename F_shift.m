function [F_f] = F_shift(csv_data,t_i_end,t_f_end,t_i_start,t_f_start,n)

F_n = csv_data.force_plates(n).F;

F_nprime = resample(F_n,round(2000*(t_f_end-t_f_start)),round(2000*(t_i_end-t_i_start)),10,40);
m = round(2000*(t_f_end-t_f_start))/round(2000*(t_i_end-t_i_start));

t_i_end_m = m*t_i_end;
shift = round(2000*(t_f_end-t_i_end_m));

t = 1:1:length(F_nprime);
t_resample = t + shift*ones(1,length(t));

F_nprime2 = F_nprime(find(t_resample > 0):1:length(t_resample),1:3);

% Chaning dimension to standard
z = t_f_end/t_i_end;
h = find(csv_data.force_plates(n).F(:,3) > 0);
i = h(1)*z;

F_f = zeros(length(csv_data.force_plates(n).F),3);
F_f(round(i):round(t_f_end*2000),:) = F_nprime2(round(i):round(t_f_end*2000),:);

end

