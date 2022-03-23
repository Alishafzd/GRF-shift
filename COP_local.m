function [localCOP] = COP_local(heel,toe,ankle,COP)
% Changing global coordinate system to local
% Defining coordinate system using hool, toe, and ankle positions
% Calculating local position of COP in the defined coordinate system

x = toe - heel;
x = x./cellfun(@norm, num2cell(x, 2));

temp2 = ankle - heel;
temp2 = temp2./cellfun(@norm, num2cell(temp2, 2));

y = cross(x, temp2);
y = y./cellfun(@norm, num2cell(y, 2));

z = cross(x, y);

R = cat(3, x, y, z);
%R = permute(R, [3 2 1]);
origin = toe;
%cor0 = cor0(1:length(origin),:);

R_cell = num2cell(R, [3 2]);
R_cell = arrayfun(@cell2mat, R_cell, 'UniformOutput', false);
R_cell = cellfun(@(x) reshape(x , 3, 3), R_cell, 'UniformOutput', false);

C = num2cell(COP-origin, 2);
func = @(A,B) (A*B');
temp = cellfun(func, R_cell, C, 'UniformOutput', false);
localCOP = [temp{:}]';

end

