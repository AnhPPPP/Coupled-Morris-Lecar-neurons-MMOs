%coupled Morris_Lecar
function dydt = CML_ode(~,y,gsyn)
dydt = zeros(4,1);
v1=y(1); w1=y(2); v2=y(3); w2=y(4);

c1=8; c2=100;
vca=120; vk=-84; vl=-60; vsyn=30;
beta=0.5; theta=-20; sigma=10; phi1=0.01; phi2=0.001;
i1=0; i2=60;
gca=4; gk=8; gl=2;
k1=-1.2; k2=18; k3=12; k4=17.4;

%Steady state functions
alpha2=1./(1+exp(-(v2-theta)/sigma));
s2=alpha2./(alpha2+beta);

minf1=0.5.*(1+tanh((v1-k1)/k2));
minf2=0.5.*(1+tanh((v2-k1)/k2));

winf1=0.5.*(1+tanh((v1-k3)/k4));
winf2=0.5.*(1+tanh((v2-k3)/k4));

tauw1=1./cosh((v1-k3)/(2*k4));
tauw2=1./cosh((v2-k3)/(2*k4));

%Define currents
ica1=gca*minf1.*(v1-vca);
ik1=gk*w1.*(v1-vk);
il1=gl*(v1-vl);
isyn1=gsyn*s2.*(v1-vsyn);

ica2=gca*minf2.*(v2-vca);
ik2=gk*w2.*(v2-vk);
il2=gl*(v2-vl);

%Equation
dydt(1)=(i1-ica1-ik1-il1-isyn1)/c1;
dydt(2)=phi1*(winf1-w1)./tauw1;
dydt(3)=(i2-ica2-ik2-il2)/c2;
dydt(4)=phi2*(winf2-w2)./tauw2;

end