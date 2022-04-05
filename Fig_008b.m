function [X_h,d_i] = Fig_008b(target_frequency)
% Fig_008b
%
% Use adjoint-derived-gradient-based optimization to maximize the growth
% rate at a target frequency by varying the heater position and iris
% diameter at downstream boundary, d_o.

%% Set the target frequency if not supplied
if nargin == 0
    target_frequency = 3.2;
end

%% Set the dimensional parameters
param_dim = fun_param_dim('Rijke');

% Increase the heater power
param_dim.Q = 500;

%% Set the numerical scheme and the max number of iterations
scheme.N     = 101;
scheme.itmax = 10;

%% Set type of adjoint, discretization, and solution method
type_of_adjoint = 'DA';
discretization  = 'FEW';
lin_or_nonlin   = 'linear'; % set to linear for robustness

%% Create the data for the contours if it doesn't already exist

if isempty(dir('Fig_008b.mat'))
    % Set the grid
    X_h = 0.0:0.01:1.0;
    d_i = 0.4:0.01:1.0;
    [XX_h,dd_i] = meshgrid(X_h,d_i);
    
    % Preallocate variables
    s      = zeros(length(d_i),length(X_h));
    ds_X_h = s;
    ds_d_i = s;

    for nn = 1:length(X_h)
        for mm = 1:length(d_i)
            param_dim.X_w = XX_h(mm,nn);
            param_dim.X_h = XX_h(mm,nn);
            % Calculate the reflection coefficient at this value of d_o
            [Rd,dRd_d_o] = fun_iris(dd_i(mm,nn));
            param_dim.Rd  = Rd;
            % Calculate the nondimensional parameters
            param = fun_nondim(param_dim);
            % Calculate the starting guess for s
            scheme.s0    = fun_set_s0(param);
            % Converge to the actual value of s
            [emode,~,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
            % Store s and ds/d(X_f) and ds/d(d_o)
            s(mm,nn) = emode.s;
            ds_X_h(mm,nn) = ds_ext.X_w + ds_ext.X_h;
            ds_d_i(mm,nn) = ds_ext.Rd * dRd_d_o;
        end
    end
    
    % Save to file
    save('Fig_008b.mat','X_h','d_i','XX_h','dd_i','s','ds_X_h','ds_d_i')
    
else
    % Load from file
    load('Fig_008b.mat')
end

%% Plot the fields and the gradients
figure(1); clf;
hold on
% plot the growth rate contours
[~,han] = contourf(XX_h,dd_i,real(s),50);
set(han,'LineColor','none'); % Switch off lines
% Plot the adjoint-based gradients
% han_q = quiver(XX_f,dd_o,real(ds_X_f),real(ds_d_o),2);
han_q = quiver(XX_h(1:2:end,1:2:end),dd_i(1:2:end,1:2:end),real(ds_X_h(1:2:end,1:2:end)),real(ds_d_i(1:2:end,1:2:end)),2);
set(han_q,'Color','w')
xlabel('Heater position, ($X_h,X_w$)','Interpreter','Latex')
ylabel('Iris diameter, $d_i$','Interpreter','Latex')
CLim = get(gca,'CLim');

% Plot the frequency contours in black
[c,han] = contour(XX_h,dd_i,imag(s),1.0:0.2:4.0);
set(han,'LineColor',[0 0 0]);
clabel(c,han);

% Format the picture
set(gca,'CLim',CLim);
axis equal; %colormap(viridis)
axis([X_h(1),X_h(end),d_i(1),d_i(end)])

set(gca,'FontSize',24)

han_cbar = colorbar;
han_cbar.Label.String = '$s_r$';
han_cbar.Label.Interpreter = 'Latex';
han_cbar.Label.FontSize = 18;

[colmap] = figfun_colmap();


% %% Perform optimization (Exterior point method with a penalty function)
% % Using a Quasi-Newton method (BFGS)
% % (from tutorial by J-D Muller, on www.fluids.ac.uk)
% 
% % Change to nonlinear for speed
% lin_or_nonlin = 'nonlin';
% 
% % Set starting point (reuse old variable names)
% X_h = 0.5;
% d_i = 0.5;
% plot(X_h,d_i,'ko','MarkerFaceColor',colmap(5,:));
% 
% % Set target value for s_i (the non-dimensional frequency)
% s_i_t = target_frequency;
% 
% % Set starting multiplier for penalty function
% r = 10;
% % Set tolerance for iteration
% g_tol = 1e-3;
% % Set eps for BFGS method
% eps_BFGS = 1e-15;
% 
% % Create col vector containing the two parameters
% x = [X_h ; d_i];
% 
% g_norm = 2*g_tol; % dummy g_norm
% nIter = 0;
% 
% while g_norm > g_tol % Iterate until close enough to minimum
%     nIter = nIter + 1;
%     % Set heater position
%     param_dim.X_w = x(1);
%     param_dim.X_h = x(1);
%     % calculate reflection coefficient from iris diameter
%     [Rd,dRd_d_o] = fun_iris(x(2));
%     param_dim.Rd  = Rd;
%     % Calculate the nondimensional parameters
%     param = fun_nondim(param_dim);
%     % Calculate the starting guess for s
%     scheme.s0    = fun_set_s0(param);
%     % Calculate the actual value for s and the gradients w.r.t external parameters
%     [emode,~,ds_ext] = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
%     % Calculate the cost function augmented by the penalty function (as a function to be minimized)
%     L = - real(emode.s) + r * (imag(emode.s) - s_i_t)^2;
%     % Calculate the gradients of s w.r.t. X_f and d_o
%     ds_X_h = ds_ext.X_w + ds_ext.X_h;
%     ds_d_i = ds_ext.Rd * dRd_d_o;
%     % Augment the gradients with the gradient of the penalty function
%     dL_dX_f = - real(ds_X_h) + r * 2 * (imag(emode.s) - s_i_t) * imag(ds_X_h);
%     dL_dd_o = - real(ds_d_i) + r * 2 * (imag(emode.s) - s_i_t) * imag(ds_d_i);
%     % Calculate the gradient as a column vector
%     g = [dL_dX_f ; dL_dd_o];
%     g_norm = norm(g) ;
%     
%     % New inverse Hessian estimate
%     if ( nIter ==1 )
%         H=eye(size(x,1));
%     else
%         % BFGS update
%         Dx=x-xpr;
%         Dg=g-gpr;
%         HDg=H*Dg;
%         if ( (Dx'*Dg) < eps_BFGS || ( HDg'*Dg) < eps_BFGS || (HDg'*HDg) < eps_BFGS )
%             H=eye(size(x,1));
%         else
%             u = Dx / (Dx'*Dg) - HDg / ( HDg' * Dg );
%             H = H ...
%                 + ( Dx * Dx' ) / ( Dx' * Dg ) ...
%                 - ( HDg * HDg' ) / ( Dg' * HDg ) ...
%                 + ( HDg' * Dg ) * (u * u');
%         end
%     end
%     
%     if ( abs( det(H) ) < 1.e-10 )
%         fprintf('singular H\n') ;
%     end
%     
%     % Set search direction, note H is an approx to the inverse Hessian.
%     p = -H*g ;
%     % Set step length
%     step = 0.01;
%     
%     % Remember current iterate before updating.
%     xpr = x ;
%     gpr = g ;
%     
%     % Update x with Armijo line search
%     [~, x, ~] = armijo_line_search(x, L, g, p, step, discretization,type_of_adjoint,lin_or_nonlin,param_dim,scheme,r,s_i_t);        % Line search
%     
%     % Plot current point
%     han = plot(x(1),x(2),'ko','MarkerFaceColor','k');
%     plot([xpr(1),x(1)],[xpr(2),x(2)],'k-')
%     drawnow
%     set(han,'Marker','.','MarkerSize',10)
%     
%     % Increase r, the multiplier for the constraint, up to some maximum
%     if fix(nIter/2) == nIter/2
%         r = min(r*1.1,100);
%     end
%     
% end
% 
% X_h = x(1);
% d_i = x(2);
% 
% plot(X_h,d_i,'ko','MarkerFaceColor',colmap(1,:));
% 
% print('-depsc2','figures/Fig_008b.eps')
% 
% end
% 
% function [Rd,dRd] = fun_iris(d_o)
% % Determine the reflection coefficient, Rd, as a function of the exit
% % orifice diameter, d_o. Also return the gradient d(Rd)/d(d_o).
% %
% % These values are based on experiments performed by Nick Jamieson for his
% % thesis (2018) with the multi-microphone method implemented by Hans Yu.
% % The results were not sufficiently rigorous to be published but serve as a
% % useful model for this optimization algorithm.
% 
% Rd  = + 0.97      * cos(pi*d_o) + 0.80 * 1i      * sin(pi*d_o);
% dRd = - 0.97 * pi * sin(pi*d_o) + 0.80 * 1i * pi * cos(pi*d_o);
% 
% end
% 
% function [step, x2, L2] = armijo_line_search(x, f, g, p, step, discretization,type_of_adjoint,lin_or_nonlin,param_dim,scheme,r,s_i_t)
% % Armijo line search
% % (from tutorial by J-D Muller, on www.fluids.ac.uk)
% %
% % x = col vector containing starting point
% % f = cost fu1nction at the starting point
% % g = gradient of cost function w.r.t. parameters at the starting point
% % p = search direction from the starting point
% % step = starting guess for step size
% %
% % The first Wolfe condition is not used because the search direction
% % is always parallel to the gradient direction in this implementation.
% 
% % Set parameter for second Wolfe condition
% eta1 = 0.1 ;
% % Set parameter for third Wolfe condition
% eta2 = 0.1 ;
% % Set the enlargement and reduction parameters
% C_enlrg = 3.0 ;
% c_reduce = 0.3 ;
% 
% % Calculate the dot product of p and g
% pg = p'*g ;
% 
% %% Calculate the function value at the predicted step
% x2 = x + step*p ;
% % Set heater position
% param_dim.X_w = x2(1);
% param_dim.X_h = x2(1);
% % calculate reflection coefficient from iris diameter
% [Rd,~] = fun_iris(x2(2));
% param_dim.Rd  = Rd;
% % Calculate the nondimensional parameters
% param = fun_nondim(param_dim);
% % Calculate the starting guess for s
% scheme.s0    = fun_set_s0(param);
% % Calculate the eigenmode
% emode = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
% % Calculate the function value
% L2 = - real(emode.s) + r * (imag(emode.s) - s_i_t)^2;
% % Calculate the ratio D
% D = ( L2-f )/step/pg ;
% 
% %% Change the step size if necessary, based on D
% sMin = 0. ;
% while D > 1-eta2 % Third Wolfe cond: enlarge step.
%     disp(['Enlarging step to ',num2str(step,12)])
%     step = min( C_enlrg*step, 0.5*step/(1-D) ) ;
%     sMin = step ;
%     x2 = x + step*p ;
%     param_dim.X_w = x2(1);
%     param_dim.X_h = x2(1);
%     [Rd,~] = fun_iris(x2(2));
%     param_dim.Rd  = Rd;
%     param = fun_nondim(param_dim);
%     scheme.s0    = fun_set_s0(param);
%     emode = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
%     L2 = - real(emode.s) + r * (imag(emode.s) - s_i_t)^2;
%     D = ( L2-f )/step/pg ;
% end
% 
% while D < eta1 && step > sMin % Second Wolfe cond: reduce step.
%     step = max( sMin + c_reduce*(step-sMin), 0.5*step/(1-D) ) ;
%     disp(['Reducing step to ',num2str(step,12)])
%     x2 = x + step*p ;
%     param_dim.X_w = x2(1);
%     param_dim.X_h = x2(1);
%     [Rd,~] = fun_iris(x2(2));
%     param_dim.Rd  = Rd;
%     param = fun_nondim(param_dim);
%     scheme.s0    = fun_set_s0(param);
%     emode = fun_Helm(discretization,type_of_adjoint,lin_or_nonlin,param,scheme);
%     L2 = - real(emode.s) + r * (imag(emode.s) - s_i_t)^2;
%     D = ( L2-f )/step/pg ;
% end

end

