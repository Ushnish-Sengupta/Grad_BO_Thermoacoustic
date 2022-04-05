function [s, pr, pl] = fun_eig_nearest(Abc,Cbc,s0)
% fun_eig_nearest
% 
% Calculate the eigenmodes of Ap = (s^2)Cp and select the one with s closest to s0
%
% INPUTS
% A  matrix representing acoustics and heat release
% C  matrix representing 1/gamma
% s0 target eigenvalue
%
% OUTPUTS
% pr right eigenvector
% s  eigenvalue
% pl left eigenvector

% Calculate the eigenvalues (s^2) and the left and right eigenvectors
[dd,ee,aa] = eig(Abc,Cbc); ee = diag(ee);

% Extract the eigenvalue closest to s^2
[ind1,d1] = dsearchn( sqrt(ee),s0); 
[ind2,d2] = dsearchn(-sqrt(ee),s0);
if abs(d1) > abs(d2)
    s = -sqrt(ee(ind2)); pr = dd(:,ind2); pl = aa(:,ind2);
else
    s = +sqrt(ee(ind1)); pr = dd(:,ind1); pl = aa(:,ind1);
end

end
