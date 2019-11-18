
function FORCING = initialize_forcing(i)
    
load('CG3_CCLM_forcing_90_101.mat')

%i = 22841;
p = 1;
% 22842: 24.4550�C --> CanopyFluxesSum: energy conservation error (3): 4.3089
% 21522: 20.9400�C --> CanopyFluxesSum: energy conservation error (3): 4.3771
% 21630: 10.4162�C --> CanopyFluxesSum: energy conservation error (3): 8.0616

% FORCING.i.Sin(p)=FORCING.data.Sin(995);
% FORCING.i.Lin(p)=FORCING.data.Lin(995);
% FORCING.i.rainfall(p)=FORCING.data.rainfall(995);
% FORCING.i.snowfall(p)=FORCING.data.snowfall(995);
% FORCING.i.q(p)=FORCING.data.q(995);
% FORCING.i.Tair(p)=FORCING.data.Tair(995);
% FORCING.i.wind(p)=FORCING.data.wind(995);
% FORCING.i.p(p)=FORCING.data.p(995);

% FORCING.i.Sin(p)=FORCING.data.Sin(22192);
% FORCING.i.Lin(p)=FORCING.data.Lin(22192);
% FORCING.i.rainfall(p)=FORCING.data.rainfall(22192);
% FORCING.i.snowfall(p)=FORCING.data.snowfall(22192);
% FORCING.i.q(p)=FORCING.data.q(22192);
% FORCING.i.Tair(p)=FORCING.data.Tair(22192);
% FORCING.i.wind(p)=FORCING.data.wind(22192);
% FORCING.i.p(p)=FORCING.data.p(22192);

% FORCING.i.Sin(p)=FORCING.data.Sin(21630);
% FORCING.i.Lin(p)=FORCING.data.Lin(21630);
% FORCING.i.rainfall(p)=FORCING.data.rainfall(21630);
% FORCING.i.snowfall(p)=FORCING.data.snowfall(21630);
% FORCING.i.q(p)=FORCING.data.q(21630);
% FORCING.i.Tair(p)=FORCING.data.Tair(21630);
% FORCING.i.wind(p)=FORCING.data.wind(21630);
% FORCING.i.p(p)=FORCING.data.p(21630);

% % % percipitation: [3.3831766e-05], temp: 21.1136
% % FORCING.i.Sin(p)=FORCING.data.Sin(22860); %(21522);
% % FORCING.i.Lin(p)=FORCING.data.Lin(22860); %(21522);
% % FORCING.i.rainfall(p)=FORCING.data.rainfall(22860); %(21522);
% % FORCING.i.snowfall(p)=FORCING.data.snowfall(22860); %(21522);
% % FORCING.i.q(p)=FORCING.data.q(22860); %(21522);
% % FORCING.i.Tair(p)=FORCING.data.Tair(22860); %(21522);
% % FORCING.i.wind(p)=FORCING.data.wind(22860); %(21522);
% % FORCING.i.p(p)=FORCING.data.p(22860); %(21522);

FORCING.i =[];
FORCING.i.Sin(p)=double(FORCING.data.Sin(i)); %(21522);
FORCING.i.Lin(p)=double(FORCING.data.Lin(i)); %(21522);
FORCING.i.rainfall(p)= double(FORCING.data.rainfall(i)); %(21522);
FORCING.i.snowfall(p)= double(FORCING.data.snowfall(i)); %(21522);
FORCING.i.q(p)=double(FORCING.data.q(i)); %(21522);
FORCING.i.Tair(p)= double(FORCING.data.Tair(i)); %(21522);
FORCING.i.wind(p)= double(FORCING.data.wind(i)); %(21522);
FORCING.i.p(p)=double(FORCING.data.p(i)); %(21522);
FORCING.i.t_span(p)=double(FORCING.data.t_span(i)); %(21522);
end