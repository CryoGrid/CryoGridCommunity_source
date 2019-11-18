function [vegetation] = filter_exposedvegp(vegetation)
%this is a dummy function since we only have trees as patch so far!!!!! :)
f = vegetation.mlcanopyinst.f;
vegetation.mlcanopyinst.p=f;
end