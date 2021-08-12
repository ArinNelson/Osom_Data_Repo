function z = roms_depth_zeta(z0,zeta,h,Vtransform)

	% Depths at this time step
    z = NaN(L,M,N);
    %z_w = NaN(L,M,N+1);
    switch Vtransform
      case 1
        for j=1:N
          z(:,:,j) = z0(:,:,j) + zeta.*(1 + z0(:,:,j)./h);  
        end
        %for j=1:(N+1)
        %  z_w(:,:,j) = z0_w(:,:,j) + zeta.*(1 + z0_w(:,:,j)./h);    
        %end
      case 2
        for j=1:N
          z(:,:,j) = zeta + (1 + zeta./h).*z0(:,:,j);
        end  
        %for j=1:(N+1)
        %  z_w(:,:,j) = zeta + (1 + zeta./h).*z0_w(:,:,j);  
        %end
    end

end