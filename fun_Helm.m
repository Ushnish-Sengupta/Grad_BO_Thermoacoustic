function [emode,ds_int,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme,TT_int,TT_ext)
% fun_Helm
%
% Calculate the eigenmodes, internal sensitivities, and external
% sensitivities with a thermoacoustic Helmholtz solver. 
%
% INPUTS
%
% discretization:
% FDS       finite difference of the strong form
% FEW       finite element of the weak form
% FDW       finite difference of the weak form
% SBP       finite difference of the strong form by Summation by Parts with the Simultaneous Approximation Term
% 
% type_of_adjoint:
% DA        discrete adjoint
% CA        continuous adjoint
%
% lin_or_nonlin:
% linear    perform the active iteration method on a series of linear eigenvalue problems
% nonlin    perform Newton iteration directly on the nonlinear eigenvalue problem
%
% param: the non-dimensional parameters
%
% scheme: the numerical scheme
%
% TT_int: list of the internal variables to be Taylor tested, e.g. {'t'}
% TT_ext: list of the external variables to be Taylor tested, e.g. {'x_f'}
%
% OUTPUTS
% 
% emode: eigenvalue, right eigenvector, left eigenvector, receptivities
%
% ds_int: gradients of the eigenvalue w.r.t. each internal parameter
% 
% ds_ext: gradients of the eigenvalue w.r.t. each external parameter

%% Read in the number of Finite elements or Chebyshev points -1
N = scheme.N;

%% Read in the point positions and matrices corresponding to the different discretizations
switch discretization
    case {'FDS','FDW'}
        [pos,mat] = mat_FD(param,N);
    case 'FEW'
        [pos,mat] = mat_FE(param,N);
    case 'SBP'
        [pos,mat] = mat_SBP(param,N);
end

%% Iterate to find the left eigenvector, pl, eigenvalue, s, and right eigenvector, pr
% Set initial s
s = scheme.s0;
switch lin_or_nonlin
    case 'linear'
        % Pre-allocate for storage of eigenmodes and weighting coefficients
        sv   = zeros(1  ,scheme.itmax);
        plv  = zeros(N+1,scheme.itmax);
        prv  = zeros(N+1,scheme.itmax);
        xiv  = zeros(1  ,scheme.itmax);
        for it = 1:scheme.itmax
            % Generate the A, C, and dAds matrices
            switch type_of_adjoint
                case 'DA'
                    switch discretization
                        case 'FEW'
                            [A,C,dAds] = mat_AC_FEW_DA(mat,param,N,s);
                        case 'FDS'
                            [A,C,dAds] = mat_AC_FDS_DA(mat,param,N,s);
                        case 'FDW'
                            [A,C,dAds] = mat_AC_FDW_DA(mat,param,N,s);
                        case 'SBP'
                            [A,C,dAds] = mat_AC_SBP_DA(mat,param,N,s);
                    end
                case 'CA'
                    switch discretization
                        case 'FEW'
                            [A,C,dAds] = mat_AC_FEW_CA(mat,param,N,s);
                        case 'FDS'
                            [A,C,dAds] = mat_AC_FDS_CA(mat,param,N,s);
                        case 'FDW'
                            [A,C,dAds] = mat_AC_FDW_CA(mat,param,N,s);
                        case 'SBP'
                            [A,C,dAds] = mat_AC_SBP_CA(mat,param,N,s);
                    end
            end
            % Find the eigenvalue and right and left eigenvectors
            [s, pr, pl] = fun_eig_nearest(A,C,s);
            % Calculate and store weighting coefficient, xi(it)
            xiv(it) = (pl' * dAds * pr) / (pl' * C * pr) / 2 / s;
            % Store eigenmode, (s, pr1, and pl1)
            sv(it) = s; prv(:,it) = pr; plv(:,it) = pl;
            % Display the dimensionless eigenvalue
            disp(['fun_Helm_',discretization,'_',type_of_adjoint,'_',lin_or_nonlin,':s= ',num2str(s,'%+.16f'),' (nondim)'])
        end
    case {'nonlin','nonlinear'}
        % Set tolerance, and dummy dels
        tol = 1e-10; dels = 2*tol;
        while abs(dels) > tol
            % Generate the A, C, and dAds matrices
            switch type_of_adjoint
                case 'DA'
                    switch discretization
                        case 'FEW'
                            [A,C,dAds] = mat_AC_FEW_DA(mat,param,N,s);
                        case 'FDS'
                            [A,C,dAds] = mat_AC_FDS_DA(mat,param,N,s);
                        case 'FDW'
                            [A,C,dAds] = mat_AC_FDW_DA(mat,param,N,s);
                        case 'SBP'
                            [A,C,dAds] = mat_AC_SBP_DA(mat,param,N,s);
                    end
                case 'CA'
                    switch discretization
                        case 'FEW'
                            [A,C,dAds] = mat_AC_FEW_CA(mat,param,N,s);
                        case 'FDS'
                            [A,C,dAds] = mat_AC_FDS_CA(mat,param,N,s);
                        case 'FDW'
                            [A,C,dAds] = mat_AC_FDW_CA(mat,param,N,s);
                        case 'SBP'
                            [A,C,dAds] = mat_AC_SBP_CA(mat,param,N,s);
                    end
            end
            % Evaluate G = A - s^2 * C and dGds
            G = A - s^2*C; dGds = dAds - 2*s*C;
            % evaluate delta s with Jacobi's formula: dels = - 1/trace(G\dGds) = - 1/trace(inv(G)*dGds);
            [QQ,RR] = qr(G); dels = - 1/trace(RR\(QQ'*dGds));
            % Update s
            s = s + dels;
            % Display the dimensionless eigenvalue
            disp(['fun_Helm_',discretization,'_',type_of_adjoint,'_',lin_or_nonlin,':s= ',num2str(s,'%+.16f'),' (nondim)'])
        end
        % Find the corresponding right and left p eigenvectors
        pr = null(G); pl = null(G');
        if isempty(pr)
            disp('G is not sufficiently converged. Decrease ''tol'' or use ''linear'''); beep
            return
        end
end

%% Calculate error in Robin boundary conditions
show_boundary_error = 'off'; % 'on' | 'off
switch show_boundary_error
    case 'on'
        if abs(param.Rd) ~= 1
            switch type_of_adjoint
                case 'DA'
                    switch discretization
                        case 'FEW'
                            error_d = param.kds * s * pr(1) - mat.D01(1,:)*pr;
                        case {'FDS','FDW'}
                            error_d = param.kds * s * pr(1) - mat.D(1,:)*pr;
                        case 'SBP'
                            error_d = param.kds * s * pr(1) - mat.D1(1,:)*pr;
                    end
                case 'CA'
                    switch discretization
                        case 'FEW'
                            error_d = conj(param.kds) * s * pr(1) - mat.D01(1,:)*pr;
                        case {'FDS','FDW'}
                            error_d = conj(param.kds) * s * pr(1) - mat.D(1,:)*pr;
                        case 'SBP'
                            error_d = conj(param.kds) * s * pr(1) - mat.D1(1,:)*pr;
                    end
            end
            disp(['Error at downstream boundary = ',num2str(norm(error_d,16))])
        end
        
        if abs(param.Ru) ~= 1
            switch type_of_adjoint
                case 'DA'
                    switch discretization
                        case 'FEW'
                            error_u = param.kus * s * pr(N+1) + mat.D01(N,:)*pr;
                        case {'FDS','FDW'}
                            error_u = param.kus * s * pr(N+1) + mat.D(N+1,:)*pr;
                        case 'SBP'
                            error_u = param.kus * s * pr(N+1) + mat.D1(N+1,:)*pr;
                    end
                case 'CA'
                    switch discretization
                        case 'FEW'
                            error_u = conj(param.kus) * s * pr(N+1) + mat.D01(N,:)*pr;
                        case {'FDS','FDW'}
                            error_u = conj(param.kus) * s * pr(N+1) + mat.D(N+1,:)*pr;
                        case 'SBP'
                            error_u = conj(param.kus) * s * pr(N+1) + mat.D1(N+1,:)*pr;
                    end
            end
            disp(['Error at upstream boundary = ',num2str(norm(error_u,16))])
        end
end

%% Read in the transformation matrix, T
switch type_of_adjoint
    case 'DA'
        switch discretization
            case 'FEW'
                [~,~,~,T] = mat_AC_FEW_DA(mat,param,N,s);
            case 'FDS'
                [~,~,~,T] = mat_AC_FDS_DA(mat,param,N,s);
            case 'FDW'
                [~,~,~,T] = mat_AC_FDW_DA(mat,param,N,s);
            case 'SBP'
                [~,~,~,T] = mat_AC_SBP_DA(mat,param,N,s);
        end
    case 'CA'
        switch discretization
            case 'FEW'
                [~,~,~,T] = mat_AC_FEW_CA(mat,param,N,s);
            case 'FDS'
                [~,~,~,T] = mat_AC_FDS_CA(mat,param,N,s);
            case 'FDW'
                [~,~,~,T] = mat_AC_FDW_CA(mat,param,N,s);
            case 'SBP'
                [~,~,~,T] = mat_AC_SBP_CA(mat,param,N,s);
        end
end

%% Normalize the eigenvectors and output the eigenmode and receptivities to a structure
emode = wrap_emode(pos,mat,T,pl,s,pr,discretization,type_of_adjoint);

%% Return if no sensitivities have been asked for
if nargout == 1
    return
end

%% Calculate the internal sensitivities (ds/d*)
switch type_of_adjoint
    case 'DA'
        switch lin_or_nonlin
            case 'linear'
                % Pre-allocate chi and set the value at scheme.itmax
                chiv = zeros(1,scheme.itmax); chiv(scheme.itmax) = 1;
                % Set the values of chi from scheme.itmax to 1
                for iti = scheme.itmax-1:-1:1; chiv(iti) = xiv(iti+1) * chiv(iti+1); end
                % Calculate dA at scheme.s0 and initialize ds_int
                switch discretization
                    case 'FEW'
                        [~,~,~,~,dA,ds_int] = mat_AC_FEW_DA(mat,param,N,scheme.s0);
                    case 'FDS'
                        [~,~,~,~,dA,ds_int] = mat_AC_FDS_DA(mat,param,N,scheme.s0);
                    case 'FDW'
                        [~,~,~,~,dA,ds_int] = mat_AC_FDW_DA(mat,param,N,scheme.s0);
                    case 'SBP'
                        [~,~,~,~,dA,ds_int] = mat_AC_SBP_DA(mat,param,N,scheme.s0);
                end
                % Calculate the denominator that is common to all sensitivities
                ip    = (plv(:,1)' * C * prv(:,1)) * 2 * s;
                % Calculate internal base state sensitivities and feedback sensitivities
                ds_int = fun_ds_DA(ds_int,chiv(1),plv(:,1),dA,prv(:,1),ip,N);
                % Subsequent iterations:
                for it = 2:scheme.itmax
                    % Calculate dA at the previous value of s
                    switch discretization
                        case 'FEW'
                            [~,~,~,~,dA] = mat_AC_FEW_DA(mat,param,N,sv(it-1));
                        case 'FDS'
                            [~,~,~,~,dA] = mat_AC_FDS_DA(mat,param,N,sv(it-1));
                        case 'FDW'
                            [~,~,~,~,dA] = mat_AC_FDW_DA(mat,param,N,sv(it-1));
                        case 'SBP'
                            [~,~,~,~,dA] = mat_AC_SBP_DA(mat,param,N,sv(it-1));
                    end
                    % Calculate the denominator that is common to all sensitivities
                    ip    = (plv(:,it)' * C * prv(:,it)) * 2 * s;
                    % Increment internal base state sensitivities and feedback sensitivities
                    ds_int = fun_ds_DA(ds_int,chiv(it),plv(:,it),dA,prv(:,it),ip,N);
                end
            case {'nonlin','nonlinear'}
                % Calculate dAds and dA at the final value of s and initialize ds_int
                switch discretization
                    case 'FEW'
                        [~,~,dAds,~,dA,ds_int] = mat_AC_FEW_DA(mat,param,N,s);
                    case 'FDS'
                        [~,~,dAds,~,dA,ds_int] = mat_AC_FDS_DA(mat,param,N,s);
                    case 'FDW'
                        [~,~,dAds,~,dA,ds_int] = mat_AC_FDW_DA(mat,param,N,s);
                    case 'SBP'
                        [~,~,dAds,~,dA,ds_int] = mat_AC_SBP_DA(mat,param,N,s);
                end
                % Calculate the denominator that is common to all sensitivities
                dGds = dAds - 2*s*C; ip = -(pl' * dGds * pr);
                % Calculate internal base state sensitivities and feedback sensitivities
                ds_int = fun_ds_DA(ds_int,1,pl,dA,pr,ip,N);
        end
    case 'CA'
        % Calculate the adjoint p eigenvector 
        pl = pr;
        % Calculate the direct eigenvalue
        s = conj(s);
        % Generate the A and C matrices of the direct problem
        switch discretization
            case 'FEW'
                [A,C] = mat_AC_FEW_DA(mat,param,N,s);
            case 'FDS'
                [A,C] = mat_AC_FDS_DA(mat,param,N,s);
            case 'FDW'
                [A,C] = mat_AC_FDW_DA(mat,param,N,s);
            case 'SBP'
                [A,C] = mat_AC_SBP_DA(mat,param,N,s);
        end
        % Calculate the direct p eigenfunction
        [~, pr, ~] = fun_eig_nearest(A,C,s);
        % Calculate ds_int with the continuous sensitivitites
        ds_int = fun_ds_CA(mat,param,T,N,pl,s,pr,discretization);
end

%% Calculate the external base state sensitivities
switch type_of_adjoint
    case 'DA'
        switch discretization
            case {'FEW','SBP'}
                ds_ext = fun_ext_int_DA(ds_int,param,pos.x0,pos.x1,N);
            case {'FDS','FDW'}
                ds_ext = fun_ext_int_DA(ds_int,param,pos.x ,pos.x ,N);
        end
    case 'CA'
        switch discretization
            case {'FDS','FDW'}
                ds_ext = fun_ext_int_CA(ds_int,param,pos.x ,pos.x ,N,mat.M  ,mat.M  );
            case 'FEW'
                ds_ext = fun_ext_int_CA(ds_int,param,pos.x0,pos.x1,N,mat.M00,mat.M11);
            case 'SBP'
                ds_ext = fun_ext_int_CA(ds_int,param,pos.x0,pos.x1,N,mat.M00,mat.M  );
        end
end

%% Perform a Taylor Test if requested (TT will be accurate for DA only)
if nargin >= 6
    if (~isempty(TT_int) || ~isempty(TT_ext))
        if strcmp(type_of_adjoint,'DA')
            fun_TT(discretization,lin_or_nonlin,param,N,s,ds_int,ds_ext,TT_int,TT_ext,scheme)
        elseif strcmp(type_of_adjoint,'CA')
            disp('Please select "DA" to perform a Taylor Test')
        end
    end
end

end
