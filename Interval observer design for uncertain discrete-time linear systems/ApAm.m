function [A_plus,A_minus] = ApAm(A)
% APAM Expresses matrix (vector/variable) A as the difference of
% two nonnegative matrices (vectors/variables) A^{+} and A^{-}.
% [A_plus,A_minus] = APAM(A)
%
% A_plus is A^{+} = max{0,A} (retain nonnegative entries and set 0 otherwise).
% A_minus is A^{-} = A^{+} - A.
% A = A^{+} - A^{-}.
%
% Example:
% A = [1 -2; -3 4];
% [A_plus,A_minus] = APAM(A)

%    Version:              2
%    Author:               Weijie Ren
%    Contact:              weijie.ren@outlook.com
%    Initial modified:     Jul. 26, 2020
%    Last modified:        Apr. 06, 2023

A_plus = zeros(size(A,1),size(A,2));
index = find(A >= 0);
A_plus(index) = A(index);
A_minus = A_plus - A;


% The following version maybe unwise, which we previously used.

% m=size(A,1);
% n=size(A,2);
% A_s=zeros(m,n);
% for i=1:m
%     for j=1:n
%         if A(i,j)>=0
%             A_s(i,j)=A(i,j);
%         else
%             A_s(i,j)=0;
%         end
%     end
% end
% A_x=A_s-A;
end

