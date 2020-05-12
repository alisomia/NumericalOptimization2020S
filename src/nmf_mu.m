function [x,info] = nmf_mu(V, rank, opt)
% Multiplicative upates (MU) for non-negative matrix factorization (NMF).
% 
% Ref: 
%   Daniel D. Lee and H. Sebastian Seung, 
%       "Algorithms for non-negative matrix factorization," NIPS 2000. 
% 
%
%