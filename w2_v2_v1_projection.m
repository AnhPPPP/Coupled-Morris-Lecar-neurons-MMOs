%Parameters
c1=80; c2=100;
vca=120; vk=-84; vl=-60; vsyn=30;
beta=0.5; theta=-20; sigma=10; phi1=0.01; phi2=0.001;
i1=0; i2=60;
gca=4; gk=8; gl=2;
k1=-1.2; k2=18; k3=12; k4=17.4;
gsyn=4.4;

%Steady state functions
alpha1= @(v1) 1./(1+exp(-(v1-theta)/sigma));
alpha2= @(v2) 1./(1+exp(-(v2-theta)/sigma));

s1= @(v1) alpha1(v1)./(alpha1(v1)+beta);
s2= @(v2) alpha2(v2)./(alpha2(v2)+beta);

minf1= @(v1) 0.5.*(1+tanh((v1-k1)/k2));
minf2= @(v2) 0.5.*(1+tanh((v2-k1)/k2));

winf1= @(v1) 0.5.*(1+tanh((v1-k3)/k4));
winf2= @(v2) 0.5.*(1+tanh((v2-k3)/k4));

tauw1= @(v1) 1./cosh((v1-k3)/(2*k4));
tauw2= @(v2) 1./cosh((v2-k3)/(2*k4));

%Define currents
ica1= @(v1) gca*minf1(v1).*(v1-vca);
ik1= @(w1,v1) gk.*w1.*(v1-vk);
il1= @(v1) gl*(v1-vl);
isyn1= @(v2,v1) gsyn*s2(v2).*(v1-vsyn);

ica2= @(v2) gca*minf2(v2).*(v2-vca);
ik2= @(v2,w2) gk.*w2.*(v2-vk);
il2= @(v2) gl*(v2-vl);

%Draw trajectory
options=odeset('BDF','on','RelTol',1e-9,'AbsTol',1e-9);
time=[0 5000];
[T3,Y3]=ode15s(@CML_ode,time,[-40 0.01 -60 0.01],options,4.4);
[T31,Y31]=ode15s(@CML_ode,time,Y3(end,:),options,4.4);

%super slow manifold
SS1=importdata('ss_sn.dat');
SS2=importdata('ss_hb.dat');

w2=linspace(-0.1, 1, 300);
v1=linspace(-70, 80, 300);
v2=linspace(-70, 80, 300);
[w2,v2,v1]=meshgrid(w2,v2,v1);

w1_cri= @(v2,v1) (i1-ica1(v1)-il1(v1)-isyn1(v2,v1))./(gk*(v1-vk));
w2_ss= @(v2)(i2-ica2(v2)-il2(v2))./(gk*(v2-vk));

der_minf = @(v1) (1/(2*k2)).*((sech((v1-k1)/k2)).^2);
der_F1_v1 = @(w2,v2,v1) 0.*w2+(((-gca.*der_minf(v1).*(v1-vca)-gca.*minf1(v1)-gl-gsyn*s2(v2)).*(v1-vk)-(i1-ica1(v1)-il1(v1)-isyn1(v2,v1)))./(gk*((v1-vk).^2)));
ver_der_F1_v1 = @(v2,v1) ((-gca.*der_minf(v1).*(v1-vca)-gca.*minf1(v1)-gl-gsyn*s2(v2)).*(v1-vk)-(i1-ica1(v1)-il1(v1)-isyn1(v2,v1)))./(gk*((v1-vk).^2));

expo = @(v2) exp((theta-v2)./sigma);
der_s2 = @(v2) (beta.*expo(v2))./(sigma.*((beta.*expo(v2)+beta+1).^2));
der_F1_v2= @(v2,v1) (-gsyn*(v1-vsyn)./(gk*(v1-vk))).*der_s2(v2);

g1= @(v2,v1) phi1*(winf1(v1)-w1_cri(v2,v1))./tauw1(v1);
f2= @(v2,w2) (1/c2)*(i2-ica2(v2)-ik2(v2,w2)-il2(v2));
v1_null = @(w2,v2,v1) g1(v2,v1)-der_F1_v2(v2,v1).*f2(v2,w2);

%the fold surfaces
L=ver_der_F1_v1(v2,v1);
%the first equation in the desingularized system
M=v1_null(w2,v2,v1);

%folded singularity curves
A=importdata('lower_folded_singularity.txt');
B=importdata('upper_folded_singularity.txt');

%upper super slow manifold
SSL=importdata('allinfo_ssl_4p4.dat');
%w2_ssl=SSL(:,4);
%v1_ssl=SSL(:,7);
%w1_ssl=SSL(:,8);
%v2_ssl=SSL(:,9);

SLU=importdata('4p4_allinfo_sslw_bifurcation.dat');
%%
figure('Position',[100 200 900 450]);
%plot trajectory
plot3(Y31(:,4), Y31(:,3), Y31(:,1),'k', 'linewidth',3);

hold on

%plot super slow manifold
%lower super slow sheet
plot3(SLU(1:5380,4),SLU(1:5380,9),SLU(1:5380,7), 'r', 'linewidth',3);
plot3(SLU(5380:10000,4),SLU(5380:10000,9),SLU(5380:10000,7), 'r:', 'linewidth',3);
%saddle  node
plot3(SLU(5400,4),SLU(5400,9),SLU(5400,7), 'bo','MarkerFaceColor', 'b', 'MarkerSize', 20 )
%Hopf point
plot3(SLU(5380,4),SLU(5380,9),SLU(5380,7), 'ro','MarkerFaceColor', 'r', 'MarkerSize', 20 )

%upper super slow sheet
%superslow manifold upper branch
plot3(SSL(1:3675,4),SSL(1:3675,9),SSL(1:3675,7), 'r', 'linewidth',3);
plot3(SSL(3676:5000,4),SSL(3676:5000,9),SSL(3676:5000,7), 'r:', 'linewidth',3);
%Hopf bifurcation
plot3(SSL(3751,4),SSL(3751,9), SSL(3751,7), 'ro','MarkerFaceColor', 'r', 'MarkerSize', 20 )
%Saddle node
plot3(SSL(3675,4),SSL(3675,9), SSL(3675,7), 'bo','MarkerFaceColor', 'b', 'MarkerSize', 20 )

%CDH
plot3(w2_ss(9.0676), 9.0676, 6.6085, 'bd', 'MarkerFaceColor', 'b', 'MarkerSize', 20 )
plot3(w2_ss(-42.9492), -42.9492, -28.612, 'bd','MarkerFaceColor', 'b', 'MarkerSize', 20 )

%folded node
plot3(0.3592, 12.0312, 6.58289, 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 20 )

%folded singularity curves
plot3(A(:,1),A(:,2),A(:,3), 'g', 'linewidth', 3)
plot3(B(:,1),B(:,2),B(:,3), 'g', 'linewidth', 3)

%the fold surface
fimplicit3(L,'FaceColor', [0.72 0.27 1],'EdgeColor','none', 'FaceAlpha', 0.5)

%the first equation in desingularized system surface
fimplicit3(M,'FaceColor', [1.0 0.0 1.0],'EdgeColor','none')
 
%phase transition
plot3(Y31(6008,4),Y31(6008,3),Y31(6008,1), 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 20)
plot3(Y31(3560,4),Y31(3560,3),Y31(3560,1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 20)

plot3(Y31(4913,4),Y31(4913,3),Y31(4913,1), 'gs','MarkerFaceColor', 'g', 'MarkerSize', 20)
plot3(Y31(18329,4),Y31(18329,3),Y31(18329,1), 'gh', 'MarkerFaceColor', 'g', 'MarkerSize', 20)

axis([0 0.6 -60 60 -60 60]);
xlabel('$w_2$','Interpreter', 'Latex','Fontsize', 20);
ylabel('$V_2$','Interpreter', 'Latex','Fontsize',20);
zlabel('$V_1$','Interpreter', 'Latex', 'Rotation', 360, 'Fontsize',20);
text(0.1,-55,-55,'$\raisebox{1pt}{\textcircled{\raisebox{-1pt}{1}}}$','Interpreter','latex','FontSize',10,'FontWeight','bold')
text(0.015,15,15,'$\raisebox{1pt}{\textcircled{\raisebox{-1pt}{2}}}$','Interpreter','latex','FontSize',10,'FontWeight','bold')
text(0.275,33,33,'$\raisebox{1pt}{\textcircled{\raisebox{-1pt}{3}}}$','Interpreter','latex','FontSize',10,'FontWeight','bold')
text(0.375,-25,-25,'$\raisebox{1pt}{\textcircled{\raisebox{-1pt}{4}}}$','Interpreter','latex','FontSize',10,'FontWeight','bold')
