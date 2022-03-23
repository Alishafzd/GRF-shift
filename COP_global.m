function [cor] = COP_global(heel,toe,ankle,cor0)
% Changing local system to global
% Defining the coordinate system using hool, toe, and ankle positions

S = min(length(toe), length(cor0));

x = toe - heel;
x = x./cellfun(@norm, num2cell(x, 2));

temp2 = ankle - heel;
temp2 = temp2./cellfun(@norm, num2cell(temp2, 2));

y = cross(x, temp2);
y = y./cellfun(@norm, num2cell(y, 2));

z = cross(x, y);

R = cat(3, x, y, z);
R = R(1:S, :, :);
%R = permute(R, [3 2 1]);
origin = toe(1:S,:);
cor0 = cor0(1:S,:);
cor0_cell = num2cell(cor0, 2);

R_cell = num2cell(R, [3 2]);
R_cell = arrayfun(@cell2mat, R_cell, 'UniformOutput', false);
R_cell = cellfun(@(x) reshape(x , 3, 3), R_cell, 'UniformOutput', false);
R_inv = cellfun(@inv, R_cell, 'UniformOutput', false);

func = @(A, B) A*B';
C = cellfun(func, R_inv, cor0_cell, 'UniformOutput', false);
cor = [C{:}]' + origin;

end

