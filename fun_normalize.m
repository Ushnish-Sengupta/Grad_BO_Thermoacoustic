function [p_out,f] = fun_normalize(p_in,M)
% fun_normalize
%
% normalize p such that int_0^1 p^2 dx = 1
%
% INPUTS
% p_in  column vector to be normalized
% M     mass matrix
%
% OUTPUTS
% p_out normalized column vector
% f     normalization factor

% Calculate normalization factor
f = sqrt(p_in.'*M*p_in);
% Normalize such that transpose(p)*S*p = 1
p_out = p_in/f;

% Multiply by -1 if the value of p at the midpoint is negative
N = length(p_out);
if N/2 == fix(N/2)
    p_mid = p_out(N/2+1);
else
    p_mid = (p_out((N-1)/2) + p_out((N+1)/2))/2;
end
p_out = p_out * sign(p_mid);
    
end
 