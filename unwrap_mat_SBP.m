function [D1,D01,M,M00,M01,V00,V11,S,Ed,Eu,I,t,h,wr,dwrdx] = unwrap_mat_SBP(mat)
% Unwrap mat_SBP
%
% Unwrap parameters for SBP-SAT method from mat structure

D1    = mat.D1;
D01   = mat.D01;
M     = mat.M;
M00   = mat.M00;
M01   = mat.M01;
V00   = mat.V00;
V11   = mat.V11;
S     = mat.S;
Ed    = mat.Ed;
Eu    = mat.Eu;
I     = mat.I;
t     = mat.t;
h     = mat.h;
wr    = mat.wr;
dwrdx = mat.dwrdx;

end
