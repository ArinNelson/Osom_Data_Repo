function x = nanfill(x)

    ii = find(isnan(x));
    if(~(isempty(ii) || numel(ii)==numel(x)))
        jj = find(~isnan(x));
        x(ii) = interp1(jj,x(jj),ii,'linear');
    end
    
end