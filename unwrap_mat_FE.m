function [D01,M11,M00,M01,V00,t,h,wr,dwrdx] = unwrap_mat_FE(mat)
% Unwrap mat_FE
%
% Unwrap mat structure for FE method

D01    = mat.D01;
M11    = mat.M11;
M00    = mat.M00;
M01    = mat.M01;
V00    = mat.V00;
t      = mat.t;
h      = mat.h;
wr     = mat.wr;
dwrdx  = mat.dwrdx;

end
