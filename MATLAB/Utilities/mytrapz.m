function [xm,xs] = mytrapz(t,x)
%--------------------------------------------------------------------------
% mytrapz - trapezoidal method of estimating mean and st. dev. of x(t)
% by Arin Nelson
% last updated 09/30/2021
% 
% usage: [x_mean, x_std] = mytrapz(t,x)
%--------------------------------------------------------------------------

    % Ensure proper dimensions
    t  = t(:); 
    x  = x(:);
    
    % Timespan
    dt = t(end)-t(1);

    % Mean from trapezoidal method
    xm = sum( ((x(2:end)+x(1:end-1))./2) .* (t(2:end)-t(1:end-1)) )/dt;
    
    % Variance from trapezoidal method
    xs = sqrt( sum( ((((x(2:end)+x(1:end-1))./2)-xm).^2) .* (t(2:end)-t(1:end-1)) )/dt );

end