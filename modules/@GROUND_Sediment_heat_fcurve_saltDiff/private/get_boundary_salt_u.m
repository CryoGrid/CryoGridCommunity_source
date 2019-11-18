function ground = get_boundary_salt_u(ground, forcing)

if forcing.TEMP.saltConcForcing == 0
    saltFlux_ub = 0;
else
    saltFlux_ub = ground.STATVAR.saltDiff(1).*(ground.STATVAR.saltConc(1) - forcing.TEMP.saltConcForcing)./ ...
        abs(ground.STATVAR.layerThick(1)/2);%abs(ground.STATVAR.midptDepth(1) - ground.STATVAR.layerDepth(1));
end

% surfaceState = forcing.TEMP.surfaceState;
% 
% switch surfaceState
%     case 1 %subaerial
%         saltFlux_ub = 0;
%     case 0 %submarine
%         saltFlux_ub = (ground.STATVAR.saltConc(1) - ground.CONST.benthicSalt)./(ground.STATVAR.layerThick(1)/2);
%     case -1 %subglacial
%         saltFlux_ub = 0;
%     otherwise %interpolated in between the phases
%         if surfaceState < 1 && surfaceState > 0 %between subaerial and submarine
%             saltFlux_ub = (1-surfaceState)*(ground.STATVAR.saltConc(1) - ground.CONST.benthicSalt)./(ground.STATVAR.layerThick(1)/2);
%         elseif surfaceState > -1 && surfaceState < 0 %between subglacial and submarine
%             saltFlux_ub = (1+surfaceState)*(ground.STATVAR.saltConc(1) - ground.CONST.benthicSalt)./(ground.STATVAR.layerThick(1)/2);
%         else
%             warning('surface state undefined!')
%         end
% end

ground.TEMP.saltFlux_ub = saltFlux_ub;
ground.TEMP.saltConc_ub = forcing.TEMP.saltConcForcing;

% if saltFlux_ub ~= 0
%     figure(3)
%     plot(1, saltFlux_ub, 'o')
%     hold on
%     drawnow
% end