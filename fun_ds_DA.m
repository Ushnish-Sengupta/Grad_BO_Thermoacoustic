function [ds_int] = fun_ds_DA(ds_int,chi,pl,dA,pr,alp,N)
% fun_ds_DA
%
% Evaluate the internal sensitivities in the discrete adjoint framework and
% add them to the existing ds structure
%
% INPUTS
% ds_int         structure containing the internal sensitivities, to be updated 
% chi            the weighting of this update
% pl             the left pressure eigenvector
% dA             structure containing the sensitivities of the A matrix
% pr             the right pressure eigenvector
% alp            the inner product <pl * dG/ds * pr>
% N              N+1 is the number of gridpoints in FD. N is the number of elements in FE
%
% OUTPUTS
% ds_int         updated structure containing the internal sensitivities

% Base state sensitivities
ds_int.n   = ds_int.n   + chi *  pl'         *  dA.n    *      pr      / alp; % row * mat * col sensitivity to n
ds_int.t   = ds_int.t   + chi * (pl' * dA.tl)   .*   (dA.tr *  pr).'   / alp; % row .* row sensitivity to t
ds_int.h   = ds_int.h   + chi * (pl' * dA.hl)    *   (dA.hr *  pr)     / alp; % row *  num sensitivity to h
ds_int.wr  = ds_int.wr  + chi * (pl' * dA.wrl)   *   (dA.wrr * pr).'   / alp; % num *  row sensitivity to wr
ds_int.v   = ds_int.v   + chi * (pl' * dA.vl)   .*   (dA.vr *  pr).'   / alp; % row .* row sensitivity to v
ds_int.ku  = ds_int.ku  + chi *  pl(N+1)'    *  dA.ku   *      pr(N+1) / alp; % num *  num sensitivity to ku
ds_int.kd  = ds_int.kd  + chi *  pl(1)'      *  dA.kd   *      pr(1)   / alp; % num *  num sensitivity to kd
% Feedback sensitivities
ds_int.mru = ds_int.mru + chi * (pl' * dA.mrul) .*  (dA.mrur * pr).'   / alp; % row .* row sensitivity to feedback from u into mass eq.
ds_int.mrp = ds_int.mrp + chi * (pl' * dA.mrpl) .*  (dA.mrpr * pr).'   / alp; % row .* row sensitivity to feedback from p into mass eq.
ds_int.fru = ds_int.fru + chi * (pl' * dA.frul) .*  (dA.frur * pr).'   / alp; % row .* row sensitivity to feedback from u into momentum eq.
ds_int.frp = ds_int.frp + chi * (pl' * dA.frpl) .*  (dA.frpr * pr).'   / alp; % row .* row sensitivity to feedback from p into momentum eq.
ds_int.qpu = ds_int.qpu + chi * (pl' * dA.qpul) .*  (dA.qpur * pr).'   / alp; % row .* row sensitivity to feedback from u into energy eq.
ds_int.qpp = ds_int.qpp + chi * (pl' * dA.qppl) .*  (dA.qppr * pr).'   / alp; % row .* row sensitivity to feedback from p into energy eq.
        
end
