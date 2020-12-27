theta_m= 00;

for theta_o = [0.05:0.01:0.95]
    
    porosity = 1 - theta_o - theta_m;
    
    c_w = 4.2e6;
    c_i = 1.9e6;
    c_m = 2e6;
    c_o = 2.5e6;
    L_sl = 3.34e8;
    T0=273.15;
    beta_interface = 1./3;
    
    alpha = 8e-4; %(Pa^-1) this is not the same alpha as we have used, it is the conversion form 1/m to 1Pa (factor 10-4) 
    m=0.19;
    
    
    alpha = 1.11e-4;
    n=1.48;
    alpha = 4e-4;
    n=2;
%     alpha = 1.49e-4;
%     n=1.25;
    m=1-1./n;
    
    % sat_waterIce_list = [0.2; 0.5; 1];
    T=[-60:0.01:0]';
    
      n=1./(1-m);  

    
    % for i=1:3
    %     sat_waterIce = sat_waterIce_list(i)
    sat_waterIce_list =[0.1:0.1:1];
    for i=1:size(sat_waterIce_list,2)  
        
        sat_waterIce = sat_waterIce_list(i);
        % sat_waterIce=1;
        
        %matrix water pressure in unfrozen state
%         mwp0 = 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);
%         mwp0 = real(mwp0);

        
        C0 = theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* c_i;
        X = -L_sl./ T0 .* beta_interface;
        L0 = porosity.* L_sl;
        

        %This is the second parameter in dimensionless equation 
        gamma = C0 ./ alpha ./ X ./ L0 ;
        
       plot(theta_o, gamma,'+')
       hold on
    end
    
end

porosity_max = 0.95;
gamma_max =  ((1-porosity_max).* c_m ) ./ alpha ./ X ./ (porosity_max.* L_sl);
porosity_min = 0.05;
gamma_min =  ((1-porosity_min).* c_o  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);

(gamma_max-gamma_min)./1000
