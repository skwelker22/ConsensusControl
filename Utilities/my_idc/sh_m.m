function B=sh_m(A,s);
% function B=sh_m(A,s);
% shift a matrix A by s (B=A+Is)

B=A+eye(size(A))*s;

