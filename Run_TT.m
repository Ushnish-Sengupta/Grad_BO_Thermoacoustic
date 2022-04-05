function [] = Run_TT(discretization,type_of_adjoint)
% Perform a Taylor test on the A matrix

%% Set the dimensional parameters
param_dim = fun_param_dim;

%% Calculate the reference scales and the nondimensional parameters
[param,~] = fun_nondim(param_dim);

%% Set the numerical scheme
N     = 11;
s0    = fun_set_s0(param);

%% Read in the point positions and matrices corresponding to the different discretizations
switch discretization
    case 'FEW'
        [~,mat] = mat_FE(param,N);
    case {'FDS','FDW'}
        [~,mat] = mat_FD(param,N);
    case 'SBP'
        [~,mat] = mat_SBP(param,N);        
end

%% Read in the A, C, and dAds matrices at s0
switch type_of_adjoint
    case 'DA'
        switch discretization
            case 'FEW'
                [A0,~,dAds] = mat_AC_FEW_DA(mat,param,N,s0);
            case 'FDS'
                [A0,~,dAds] = mat_AC_FDS_DA(mat,param,N,s0);
            case 'FDW'
                [A0,~,dAds] = mat_AC_FDW_DA(mat,param,N,s0);
            case 'SBP'
                [A0,~,dAds] = mat_AC_SBP_DA(mat,param,N,s0);
        end
    case 'CA'
        switch discretization
            case 'FEW'
                [A0,~,dAds] = mat_AC_FEW_CA(mat,param,N,s0);
            case 'FDS'
                [A0,~,dAds] = mat_AC_FDS_CA(mat,param,N,s0);
            case 'FDW'
                [A0,~,dAds] = mat_AC_FDW_CA(mat,param,N,s0);
            case 'SBP'
                [A0,~,dAds] = mat_AC_SBP_CA(mat,param,N,s0);
        end
end

% Set increment size (small)
epsilon = 1e-3;
ddA = zeros(1,4);

% Calculate dA at new values of s
for nn = 1:4
    switch type_of_adjoint
        case 'DA'
            switch discretization
                case 'FEW'
                    A = mat_AC_FEW_DA(mat,param,N,s0+nn*epsilon);
                case 'FDS'
                    A = mat_AC_FDS_DA(mat,param,N,s0+nn*epsilon);
                case 'FDW'
                    A = mat_AC_FDW_DA(mat,param,N,s0+nn*epsilon);
                case 'SBP'
                    A = mat_AC_SBP_DA(mat,param,N,s0+nn*epsilon);
            end
        case 'CA'
            switch discretization
                case 'FEW'
                    A = mat_AC_FEW_CA(mat,param,N,s0+nn*epsilon);
                case 'FDS'
                    A = mat_AC_FDS_CA(mat,param,N,s0+nn*epsilon);
                case 'FDW'
                    A = mat_AC_FDW_CA(mat,param,N,s0+nn*epsilon);
                case 'SBP'
                    A = mat_AC_SBP_CA(mat,param,N,s0+nn*epsilon);
            end
    end
    % Work out dAds by Finite difference
    dA_FD = (A - A0);
    dA_AD = dAds * (nn*epsilon);
    % Find the difference between dA_FD and dA_AD via adjoints
    ddA(nn) = norm(dA_FD - dA_AD);
end

figure(10)
plot([0 1 4 9 16],[0,ddA],'k.-')

end

    
