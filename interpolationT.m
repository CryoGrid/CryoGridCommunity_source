


stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(1) = interp1(grid.MIDPOINTS,stratigraphy2.STATVAR.T,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(2) = interp1(grid.MIDPOINTS,stratigraphy2.STATVAR.T,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.t_soisno(3) = interp1(grid.MIDPOINTS,stratigraphy2.STATVAR.T,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');

stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_vol(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.water,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(1) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(1),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(2) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(2),'nearest');
stratigraphy1.STATVAR.vegetation.soilvar.h2osoi_ice(3) = interp1(stratigraphy2.STATVAR.midPoint,stratigraphy2.STATVAR.ice,stratigraphy1.STATVAR.vegetation.soilvar.zi(3),'nearest');