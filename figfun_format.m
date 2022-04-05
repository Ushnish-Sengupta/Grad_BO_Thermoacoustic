function [] = figfun_format(xlab,ylab,titl,type)
% figfun_format
%
% Format the axis
%
% INPUT
% xlab      string for the horizontal axis
% ylab      string for the vertical axis
% titl      string for the title
% type      The type of format being requested

switch type
    case 'type_1'
        FontSize = 18; box on; grid on;
    otherwise
        disp(['type ',type,' not recognized in figfun_format'])
end

set(gca,'FontSize',FontSize,'FontName','Times')

if xlab
    xlabel(xlab,'Interpreter','Latex','FontSize',FontSize)
end

if ylab
    ylabel(ylab,'Interpreter','Latex','FontSize',FontSize)
end

if titl
    title(titl,'Interpreter','Latex','FontSize',FontSize)
end

end
