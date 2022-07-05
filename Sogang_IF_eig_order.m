% ====================================================================== %
% eig_order.m
% 
% given a square matrix A, this program returns a diagonal matrix L of 
% eigenvalues and a matrix X of eigenvectors where L has been sorted in
% ascending modulus. The index 'j' is the first element of L that is greater
% than 1 (in abs value).
%
% code ... [L X j] = eig_order(A) --> No commas (why?)
%
% C. House 12/18/01
% ====================================================================== %

function [L,X,j] = eig_order(A);

[x1,L1] = eig(A);		% x1 matrix of eigenvectors -- columns

L2 = diag(L1); 			% column of eigenvalues

n = length(L2);

[L3,k] = sort(abs(L2));		% k is an index vector describing the sort proceedure.

j = min(find(abs(L3)>1));

L4 = [];
X  = [];

for m = 1:n;
	L4 = [L4; L2(k(m))];
	X  = [X x1(:,k(m))];
end; 

L = diag(L4);
