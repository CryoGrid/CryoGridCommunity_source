
function FORCING = initialize_forcing(i)
    
load('CG3_CCLM_forcing_90_101.mat')

%i = 22841;
p = 1;


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