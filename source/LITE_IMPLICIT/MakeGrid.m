function [Zn, Zs, dxn, dxs, kn, ks] = MakeGrid(Zp,dxp,kp)

N = length(Zp);
Zn = Zp-dxp/2;
Zs = Zp+dxp/2;

dxn=ones(N,1);
dxs=ones(N,1);
for i=2:N-1
    dxn(i,1) = Zp(i) - Zp(i-1);
    dxs(i,1) = Zp(i+1) - Zp(i);
end

dxs(1,1) = Zp(2) - Zp(1);
dxs(N,1) = dxp(N)/2;

dxn(1,1) = dxp(1,1)/2;
dxn(N,1) = Zp(N) - Zp(N-1);

kn=ones(N,1);
ks=ones(N,1);
for i=2:N-1
    kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
    ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
end
kn(1,1) = kp(1);
kn(N,1) = kp(N);
ks(1,1) = kp(1);
ks(N,1) = kp(N);

end
