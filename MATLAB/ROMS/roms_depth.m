function z0_r = roms_depth(h,N,opts,varargin)
%=========================================================================%
% z = roms_depth(h,n,opts)
% computes depth levels given 2D bathymetry h, the number of vertical
% levels n, and the grid options in structure opts containing the following
% variables:
%   Vstretching
%   Vtransform
%   theta_s
%   theta_b
%   hc
% 
% z = roms_depth(h,s,opts,zeta)
% includes 2D free-surface zeta

% Grid options
Vstretching = opts.Vstretching;
Vtransform  = opts.Vtransform;
theta_s     = opts.theta_s;
theta_b     = opts.theta_b;
hc          = opts.hc;

% Grid dimensions
L = size(h,1);
M = size(h,2);

% The sigma levels
%k_w = 0:N;
k_r = (1:N)-0.5;
if(Vstretching<=4)
  %s_w = (k_w-N) ./ N;
  s_r = (k_r-N) ./ N;
else
  %s_w   = -(( k_w.^2 - 2*N.*k_w + k_w + N^2 - N)./(N^2 - N)) - 0.01.*( (k_w.^2 - N.*k_w) ./ (1-N) );
  s_r   = -(( k_r.^2 - 2*N.*k_r + k_r + N^2 - N)./(N^2 - N)) - 0.01.*( (k_r.^2 - N.*k_r) ./ (1-N) );
end

% Compute the vertical stretching function
switch Vstretching
  case 1
    Cs_r = (1-theta_b).*( sinh(theta_s.*s_r)./sinh(theta_s) ) + theta_b.*( (tanh(theta_s.*(s_r+0.5))./(2*tanh(0.5*theta_s))) - 0.5 ); 
    %Cs_w = (1-theta_b).*( sinh(theta_s.*s_w)./sinh(theta_s) ) + theta_b.*( (tanh(theta_s.*(s_w+0.5))./(2*tanh(0.5*theta_s))) - 0.5 ); 
  case {2,3}
    gamma = 3;  % Set internally, see wiki page on Vertical_S-coordinate
    mu_r     = 0.5.*( 1 - tanh( gamma.*( s_r+0.5 ) ) );
    %mu_w     = 0.5.*( 1 - tanh( gamma.*( s_w+0.5 ) ) );
    Cs_srf_r =  0 - log(cosh( gamma.*( abs(s_r).^theta_s ) ))./log(cosh(gamma));
    Cs_bot_r = -1 + log(cosh( gamma.*(    (s_r).^theta_b ) ))./log(cosh(gamma));
    %Cs_srf_w =  0 - log(cosh( gamma.*( abs(s_w).^theta_s ) ))./log(cosh(gamma));
    %Cs_bot_w = -1 + log(cosh( gamma.*(    (s_w).^theta_b ) ))./log(cosh(gamma));
    Cs_r     = mu_r.*Cs_bot_r + (1-mu_r).*Cs_srf_r;
    %Cs_w     = mu_w.*Cs_bot_w + (1-mu_w).*Cs_srf_w;
  case {4,5}
    if(theta_s > 0)
      Cs_0_r = (1-cosh(theta_s.*s_r)) ./ (cosh(theta_s)-1);   
      %Cs_0_w = (1-cosh(theta_s.*s_w)) ./ (cosh(theta_s)-1); 
    else
      Cs_0_r = -(s_r.^2);
      %Cs_0_w = -(s_w.^2);
    end
    Cs_r = (exp(theta_b.*Cs_0_r) - 1) ./ (1-exp(-theta_b));
    %Cs_w = (exp(theta_b.*Cs_0_w) - 1) ./ (1-exp(-theta_b));
  otherwise error(['Unknown vertical stretching function choice: ' num2str(Vstretching)]);
end

% Perform vertical transform for free-surface case
z0_r = zeros(L,M,N);
%z0_w = zeros(L,M,N+1);
switch Vtransform
    case 1
      for i=1:N
        z0_r(:,:,i) = hc*s_r(i) + (h-hc).*Cs_r(i);
      end
      %for i=1:(N+1)
      %  z0_w(:,:,i) = hc*s_w(i) + (h-hc).*Cs_w(i);
      %end
    case 2
      for i=1:N
        z0_r(:,:,i) = h.*( (hc*s_r(i) + h.*Cs_r(i)) ./ (hc + h) ); 
      end
      %for i=1:(N+1)
      %  z0_w(:,:,i) = h.*( (hc*s_w(i) + h.*Cs_w(i)) ./ (hc + h) ); 
      %end
    otherwise
    error(['Possible valid Vtransform values are 1 and 2. Given: ' num2str(Vtransform)]);
end

end