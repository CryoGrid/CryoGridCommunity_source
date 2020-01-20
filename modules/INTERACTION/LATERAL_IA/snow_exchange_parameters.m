function factor = snow_exchange_parameters(lateral)
    status = lateral.STATUS.snow;
    alt = lateral.TEMP.surfaceAltitudes;
    delta = lateral.PARA.delta;
    area = lateral.PARA.area;
    mass = lateral.TEMP.waterIce;
    in = zeros(size(alt));
    
    for i = 1:length(alt)
        I = find(alt + delta <= alt(i) & status == 1);
        if isempty(I)
            in(i) = in(i) + area(i).*mass(i);
        else
            out = area(i).*mass(i);
            in(I) =in(I) + (out.*area(I))./sum(area(I));
        end
    end
    
    factor = in./sum(mass.*area);
    
    if sum(factor) ~= 1
        warning(['Snow exchange is not mass conservingby a factor of ' num2str( sum(factor) ) '!'])
    end
end
