function L = ichol0(A)
% ICHOL0  Incomplete Cholesky factorization with zero fill-in of a
% symmetric positive definite sparse matrix A (Sec. 2.5 of the project:
% "Incomplete Choleski decomposition with zero fill-in").
%
% From Matlab R2015a+ this is simply a wrapper around the built-in
% ichol, which supports the 'nofill' type natively and is far more
% efficient than a hand-rolled loop.

opts.type = 'nofill';
opts.michol = 'off';
L = ichol(A, opts);
end
