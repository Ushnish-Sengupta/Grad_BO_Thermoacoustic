function [D,M,V,t,h,wr,dwrdx] = unwrap_mat_FD(mat)
% Unwrap mat_FD
%
% Unwrap mat structure for FD matrices

D     = mat.D;
M     = mat.M;
V     = mat.V;
t     = mat.t;
h     = mat.h;
wr    = mat.wr;
dwrdx = mat.dwrdx;

end
