function ground = initializeExcessIce(ground) 

ground.STATVAR.excessGroundIce = ground.STATVAR.waterIce ./ ground.STATVAR.layerThick > ground.STATVAR.naturalPorosity;