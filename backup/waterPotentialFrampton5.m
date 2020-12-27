
    
c_w = 4.2e6;
c_i = 1.9e6;
c_m = 2e6;
c_o = 2.5e6;
L_sl = 3.34e8;
T0=273.15;
beta_interface = 1./3;

T_min = -60;
sat_waterIce_min = 0.005;

alpha = 8e-4; %(Pa^-1) this is not the same alpha as we have used, it is the conversion form 1/m to 1Pa (factor 10-4)
m=0.19;


    alpha = 1.11e-4;
    n=1.48;
alpha = 4e-4;
n=2;
%     alpha = 1.49e-4;
%     n=1.25;
m=1-1./n;
T_min = -60;

T=[T_min:0.01:0]';


n=1./(1-m);


LUT_gamma_size = 512;
porosity_max = 0.95;
porosity_min = 0.05;

X = -L_sl./ T0 .* beta_interface;
gamma_max =  ((1-porosity_max).* c_o ) ./ alpha ./ X ./ (porosity_max.* L_sl);
log(-gamma_max)

gamma_min =  ((1-porosity_min).* c_m  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);
log(-gamma_min)
gamma_list = linspace(log(-gamma_min), log(-gamma_max), LUT_gamma_size);
gamma_list = -exp(gamma_list);

%dimensionless quantities computed from the dimensionless equation

LUT_T_prime = [];
LUT_E_prime = [];

for gamma_index=1:size(gamma_list,2)
    gamma= gamma_list(1,gamma_index);
    
    T_prime_max = alpha .* X .* T_min + ((sat_waterIce_min.^(-1./m)-1)).^(1./n);
    E_prime_min = gamma .* T_prime_max + (1 + T_prime_max.^n).^(-m);
    
    T_prime_max = T_prime_max .*1.1;
    E_prime_min = E_prime_min .*1.1;
    
    
    LUT_size = 2.^10;
    
    T_prime2 = linspace(T_prime_max, 0, LUT_size+1)';
    E_prime2 = gamma .* T_prime2 + (1 + T_prime2.^n).^(-m);
    

    for count = 1:15
        %             plot(E_prime2)
        %             hold on
        %add points
        add_points_T =[(T_prime2(1:end-1,1) + T_prime2(2:end,1))./2  (3.*T_prime2(1:end-1,1) + T_prime2(2:end,1))./4 (T_prime2(1:end-1,1) + 3.*T_prime2(2:end,1))./4];%midpoints in T
        add_points_E = gamma .* add_points_T + (1 + add_points_T.^n).^(-m);
        
        delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
        delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
        add_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (add_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
        
        add_points_T = add_points_T(:);
        add_points_E = add_points_E(:);
        mismatch = abs(add_points_T_interp(:) - add_points_T);
        [~,sort_index]=sort(mismatch);
        sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
        T_prime2 = [T_prime2; add_points_T(sort_index,1)];
        E_prime2 = [E_prime2; add_points_E(sort_index,1)];
        [T_prime2, sort_index] = sort(T_prime2,'descend');
        E_prime2=E_prime2(sort_index,1);
        
        %reduce points
        reduce_points_T = [T_prime2(2:4:end-3,1) T_prime2(3:4:end-2,1) T_prime2(4:4:end-1,1)];
        reduce_points_E = [E_prime2(2:4:end-3,1) E_prime2(3:4:end-2,1) E_prime2(4:4:end-1,1)];
        T_prime2 = T_prime2(1:4:end,1);
        E_prime2 = E_prime2(1:4:end,1);
        
        delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
        delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
        reduce_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (reduce_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
        
        reduce_points_T = reduce_points_T(:);
        reduce_points_E = reduce_points_E(:);
        mismatch = abs(reduce_points_T_interp(:) - reduce_points_T);
        [~,sort_index]=sort(mismatch);
        sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
        T_prime2 = [T_prime2; reduce_points_T(sort_index,1)];
        E_prime2 = [E_prime2; reduce_points_E(sort_index,1)];
        [T_prime2, sort_index] = sort(T_prime2,'descend');
        E_prime2 = E_prime2(sort_index,1);
        

        dT_prime2_dgamma = -(E_prime2 - (1+T_prime2.^n).^(-m)) ./ T_prime2.^-2 + n.*m .* T_prime2.^(n-2) .* (1+T_prime2.^n).^(-m-1);
        % dT_prime2_dgamma = 1./dT_prime2_dgamma;
        
    end
    
    LUT_T_prime = [LUT_T_prime T_prime2];
    LUT_E_prime = [ LUT_E_prime E_prime2];
    
end

%%
close all
        %---------------interpolate and test
        test_T = [[-50:0.01:-0.01] -0.001];
        sat_waterIce = 0.005;
        theta_o= 0.05;
        theta_m = 0;
        porosity = 1 - theta_o - theta_m;

        mwp0 = 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);

        C0 = theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* c_i;
        X = -L_sl./ T0 .* beta_interface;
        L0 = porosity.* L_sl;

        gamma = C0 ./ alpha ./ X ./ L0 ;
        
        
        mpw_test = -L_sl.*test_T./T0 .*beta_interface .* double(test_T<0) + mwp0;
        sat_water_test = double(mpw_test>0) .* (1+(alpha.*mpw_test).^n).^(-m) + double(mpw_test<=0);
        sat_ice_test = sat_waterIce - sat_water_test;
        energy_test = test_T.* (theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* (c_w .* double(test_T>=0)+ c_i.* double(test_T<0)));
        energy_test = energy_test - double(test_T<0) .* L_sl .* porosity .* (sat_waterIce - sat_water_test);
        E_prime_test = energy_test ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; 
        T_prime_test = alpha .* (X.*test_T + mwp0);

        left_gamma_index=find(gamma-gamma_list>=0, 1, 'last');
        scale_factor_gamma = (gamma-gamma_list(left_gamma_index) ) ./ (gamma_list(left_gamma_index+1) - gamma_list(left_gamma_index));

        T_prime_left = interp1(LUT_E_prime(:,left_gamma_index), LUT_T_prime(:,left_gamma_index), E_prime_test, 'linear');
        T_test = (T_prime_left ./ alpha - mwp0) ./ X;
        
%         plot(test_T, test_T)
%         hold on
%         plot(test_T, T_test)
         
        T_prime_right = interp1(LUT_E_prime(:,left_gamma_index+1), LUT_T_prime(:,left_gamma_index+1), E_prime_test, 'linear');
        T_test = (T_prime_right ./ alpha - mwp0) ./ X;
%         plot(test_T, T_test)
        
        T_prime_final = T_prime_left + scale_factor_gamma .*( T_prime_right - T_prime_left);
        T_final = (T_prime_final ./ alpha - mwp0) ./ X;
        sat_water_final = E_prime_test - gamma .* T_prime_final;
        plot(sat_water_test, sat_water_test - sat_water_final)
%         plot(test_T, T_final)
        
        figure
         plot(test_T, test_T - T_final)
        
        
        T_prime0 = alpha .*mwp0;
        T_prime_alternative = (E_prime_test - (1 + T_prime0.^n).^(-m) - T_prime0 .* (n.*m).* T_prime0.^(n-1) .* (1+T_prime0.^n).^(-m-1)) ./ (gamma -  n.*m.* T_prime0.^(n-1) .* (1+T_prime0.^n).^(-m-1));     
        T_final_alternative = (T_prime_alternative ./ alpha - mwp0) ./ X;
        sat_water_final_alternative = E_prime_test - gamma .* T_prime_alternative;
        figure
        plot(sat_water_test, sat_water_test - sat_water_final_alternative)
        
        figure
        plot(test_T, test_T - T_final_alternative)
        
        %plot((test_T .* X + mwp0).*alpha, T_prime_test)
%         plot(E_prime_test, T_prime_test)
%         hold on
        
        %
       
        
        %dT_dgamma = interp1(E_prime2, dT_prime2_dgamma, E_prime_test, 'linear');
        %dg=[dg; dT_dgamma];
%         dg=[dg; T_prime_test ./(gamma - n.*m.*T_prime_test.^(n-1) .* (1+T_prime_test.^n).^(-m-1) ) ];
%         dT=[dT; T_prime_test];
%          T_prime_test = T_prime_test + ind.* 3.7e-6 ./2 .* T_prime_test ./(gamma - n.*m.*T_prime_test.^(n-1) .* (1+T_prime_test.^n).^(-m-1) ) ;
% %         
%          T_test = (T_prime_test ./ alpha - mwp0) ./ X;
% % 
%          plot(test_T, T_test)
%          hold on
%          
%         T_acc = T_acc + T_test;
  
        
%          plot(test_T, test_T - T_test)
%          hold on
%          plot(test_T, T_acc/2)
%          hold on

% theta_o= 0.5;
% 
% for theta_m = 0;%0.95; %[0.05:0.05:0.95]
    
%     porosity = 1 - theta_o - theta_m;
% 
%         mwp0 = 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);
%         
% 
%         C0 = theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* c_i;
%         X = -L_sl./ T0 .* beta_interface;
%         L0 = porosity.* L_sl;
%         
%         %dimensionless quantities computed from the "true" quantities
% 
% 
%         gamma0 = C0 ./ alpha ./ X ./ L0 ;

% porosity_max = 0.95;
% gamma_max =  ((1-porosity_max).* c_o ) ./ alpha ./ X ./ (porosity_max.* L_sl);
% log(-gamma_max)
% porosity_min = 0.05;
% gamma_min =  ((1-porosity_min).* c_m  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);
% log(-gamma_min)


%
