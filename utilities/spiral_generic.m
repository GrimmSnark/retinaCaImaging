function A = spiral_generic(n, P)
% Makes NxN matrix filled up spirally starting with point P
r = max([P - 1, n - P]);              % Radius of the bigger matrix
M = spiral(2 * r + 1);                % Bigger matrix itself
M = permute(M,[2,1]);                 % changing start direction of the spiral
M = M(:,end:-1:1);                    % chaning spin orientation
C = r + 1 - (P - 1);                  % Top-left corner of A in M
A = M(C(1):C(1)+n(1)-1, C(2):C(2)+n(2)-1);  % Get the submatrix
[~, order] = sort(A(:));              % Get elements' order
A(order) = 1:(n(1)*n(2));                     % Fill with continous values
end