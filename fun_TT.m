function [] = fun_TT(discretization,lin_or_nonlin,param,N,s,ds_int,ds_ext,TT_int,TT_ext,scheme)
% fun_TT
%
% Perform a Taylor Test on the fields held in TT_int and TT_ext. 
% This function plots a black line and a straight grey line.
%
% If the black line lies on top of the grey line then the sensitivities 
% are correct to first order, as required. 
%
% If the black line bends like sqrt(epsilon^2) then the sensitivities are
% incorrect to first order and there is a bug in the adjoint code.
%
% If the black line bends like (epsilon^2)^{>1} then 2nd and 3rd order
% effects are influential over the range of d(parameter) being considered.
% epsilon in fun_rand_fields should be reduced. 
%
% If the line is jagged and the value on the vertical axis is small then
% increase epsilon in fun_rand_fields. If the line remains jagged then 
% the sensitivity is probably zero and the plot is showing machine errors
% 
% The Taylor Tests are performed on random perturbation vectors so if the
% line is just slightly bent, it is worth running the test again to see if
% the result changes.
%
% If the line is just slightly bent, try increasing the accuracy at which
% the eigenmodes are calculated. For nonlin this is done by setting tol to
% a lower value in fun_Helm. For linear this is done by increasing the
% number of iterations in scheme.itmax.
%
% INPUTS
% discretization discretization (FDS, FEW,, FDW, or SBP)
% lin_or_nonlin  iteration procedure (lin = active iteration; nonlin = Newton method) 
% param          structure containing the parameters
% N              N+1 = number of gridpoints for FD. N = number of elements for FE 
% s              eigenvalue
% ds_int         sensitivities of s w.r.t. internal parameters
% ds_ext         sensitivities of s w.r.t. external parameters
% TT_int         cell array of internal parameters to vary in the Taylor Test
% TT_ext         cell array of external parameters to vary in the Taylor Test
% scheme.s0      starting value of s used in the active iteration method 
%
% OUTPUTS
% none

% If TT_int = 'all' then test all the internal parameters
if strcmp(TT_int,'all')
    TT_int = {'t','h','wr','v','ku','kd'};
end

% If TT_ext = 'all' then test all the external parameters
if strcmp(TT_ext,'all')
    TT_ext = {'n','X_w','L_w','X_h','X_h','tau'};
    if param.Ru ~= +1 && param.Ru ~= -1
        TT_ext(length(TT_ext)+1) = {'Ru'};
    end
    if param.Rd ~= +1 && param.Rd ~= -1
        TT_ext(length(TT_ext)+1) = {'Rd'};
    end
end

% Label the perturbation vectors
d_int = ds_int;
d_ext = ds_ext;

% Set all perturbation vectors to zero
d_int = fun_zero_fields(d_int);
d_ext = fun_zero_fields(d_ext);

% Set chosen fields' perturbation vectors to random numbers of size arg(3)
d_int = fun_rand_fields(d_int,TT_int,1e-1);
d_ext = fun_rand_fields(d_ext,TT_ext,1e-3);

% Calcuate ds with adjoints
ds_AD = fun_add_fields(0    ,d_int,ds_int);
ds_AD = fun_add_fields(ds_AD,d_ext,ds_ext);

% Store old external parameter values and s
s_old = s;
param_old = param;

% Perform the Taylor test
ds_FD = zeros(1,4);
for nn = 1:4
    % Update external parameters
    param = fun_update_ext(param_old,d_ext,nn);
    % Update kus and kds
    [kus,kds,~,~] = fun_kukd(param);
    param.kus = kus;
    param.kds = kds;
    % Update the internal parameters
    switch discretization
        case 'FEW'
            [~,mat] = mat_FE(param,N);
            mat.V00 = mat.V00 + diag(nn*d_int.v);
        case {'FDS','FDW'}
            [~,mat] = mat_FD(param,N);
            mat.V   = mat.V   + diag(nn*d_int.v);
        case 'SBP'
            [~,mat] = mat_SBP(param,N);
            mat.V00 = mat.V00 + diag(nn*d_int.v);
    end
    mat.t   = mat.t   + nn*d_int.t;
    mat.h   = mat.h   + nn*d_int.h;
    mat.wr  = mat.wr  + nn*d_int.wr;
    param.kus = param.kus + nn*d_int.ku;
    param.kds = param.kds + nn*d_int.kd;
    % Iterate to find new s with FD
    switch lin_or_nonlin
        case {'nonlin','nonlinear'}
            tol = 1e-8; dels = 2*tol;
            while abs(dels) > tol
                switch discretization
                    % Generate the A, C, and dAds matrices
                    case 'FEW'
                        [A,C,dAds] = mat_AC_FEW_DA(mat,param,N,s);
                    case 'FDS'
                        [A,C,dAds] = mat_AC_FDS_DA(mat,param,N,s);
                    case 'FDW'
                        [A,C,dAds] = mat_AC_FDW_DA(mat,param,N,s);
                    case 'SBP'
                        [A,C,dAds] = mat_AC_SBP_DA(mat,param,N,s);
                end
                % Evaluate G = A - s^2 * C and dGds
                G = A - s^2*C; dGds = dAds - 2*s*C;
                % evaluate delta s with Jacobi's formula: dels = - 1/trace(L\dLds)
                [QQ,RR] = qr(G); dels = - 1/trace(RR\(QQ'*dGds));
                % Update s
                s = s + dels;
            end
        case 'linear'
            s = scheme.s0;
            for it = 1:scheme.itmax
                switch discretization
                    % Generate the A and C matrices
                    case 'FEW'
                        [A,C] = mat_AC_FEW_DA(mat,param,N,s);
                    case 'FDS'
                        [A,C] = mat_AC_FDS_DA(mat,param,N,s);
                    case 'FDW'
                        [A,C] = mat_AC_FDW_DA(mat,param,N,s);
                    case 'SBP'
                        [A,C] = mat_AC_SBP_DA(mat,param,N,s);
                end
                s = fun_eig_nearest(A,C,s);
            end
        otherwise
            disp('please specify linear or nonlin/nonlinear')
    end
    ds_FD(nn) = s - s_old;
end
ds_TT = ds_FD - ds_AD*[1,2,3,4];
figure(1); cla; hold on
plot([0,16],[0,abs(ds_TT(4))],'-','LineWidth',6,'Color',[0.8 0.8 0.8])
plot([0,1,4,9,16],[0,abs(ds_TT)],'k.-')
xlabel('$\epsilon^2$','Interpreter','Latex')
ylabel('$\delta s_{\mbox{adjoint}} - \delta s_{\mbox{finite difference}}$','Interpreter','Latex')
end

%% Utility functions

function d = fun_zero_fields(d)

names = fieldnames(d);
for nn = 1:length(names)
    f = getfield(d,names{nn});
    d = setfield(d,names{nn},0*f');
end

end

function d = fun_rand_fields(d,TT,epsilon)

for nn = 1:length(TT)
    f = getfield(d,TT{nn});
    d = setfield(d,TT{nn},epsilon*(rand(size(f))-0.5));
end

end

function ds_AD = fun_add_fields(ds_AD,d,ds)

names = fieldnames(d);
for nn = 1:length(names)
    ds_AD = ds_AD + getfield(ds,names{nn}) * getfield(d,names{nn});
end

end

function param = fun_update_ext(param_old,d_ext,n)

param = param_old;
names = fieldnames(param_old);
for nn = 1:length(names)
    f = getfield(param_old,names{nn});
    g = getfield(d_ext,names{nn});
    param = setfield(param,names{nn},f+g*n);
end

end

