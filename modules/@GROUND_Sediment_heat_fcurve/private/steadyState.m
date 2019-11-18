function T_0 = steadyState(ground)
%   steadyState
%
%

% get thermal conductivities
k_water = ground.CONST.k_w;
k_ice = ground.CONST.k_i;
k_organic = ground.CONST.k_o;
k_mineral = ground.CONST.k_m;

if isempty(ground.PREVIOUS)
    TForcing = ground.TEMP.T_ub;
else
    TForcing = ground.PREVIOUS.STATVAR.T(end);
end

Q = ground.PARA.heatFlux_lb;

midpoints = ground.STATVAR.upperPos - ground.STATVAR.layerThick(1)/2 - ...
    cumsum(ground.STATVAR.layerThick);
    
Tmelt = ground.STATVAR.Tmelt; %in Kelvin
Tmelt_inDegreeC = Tmelt - 273.15;%convert to °C to compare against temperatures.
    
% allocate memory
T_0 = zeros(length(midpoints),1);
% apply forcing to surface
T_0(1) = TForcing;


for i=2:size(T_0,1)
    T = T_0(i-1);

    a1=1./ground.STATVAR.porosity(i) - ground.PARA.a(i).*Tmelt(i) - ...
        ground.PARA.b(i).*Tmelt(i).^2;

    licWater = (T<Tmelt_inDegreeC(i)).*(1./(a1 + ground.PARA.a(i).*T +...
        ground.PARA.b(i).*T.^2)) + (T>=Tmelt_inDegreeC(i)).*ground.STATVAR.porosity(i);

    thermCond =(licWater.*sqrt(k_water)+(ground.STATVAR.porosity(i)-licWater).*sqrt(k_ice) + ...
        ground.STATVAR.mineral(i)*sqrt(k_mineral) + ground.STATVAR.organic(i)*sqrt(k_organic)).^2;

	T_0(i,1)=T_0(i-1,1)- 1./thermCond.*Q.*(midpoints(i,1)- midpoints(i-1,1));
end

figure(1)
plot(T_0,midpoints)
hold on
drawnow
