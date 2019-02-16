close all
clear

R=transpose(linspace(0,1,5001));
N=transpose(0.05:0.05:2);

K0=0;
KF=5e5;
Kstep=0.05;
ORDmax=20;
ORD=20;
ResLayer=500;

for d=2:1:6
    G=gwlc(R,N,d,K0,KF,Kstep,ORDmax,ORD,ResLayer);
    
    save(['data/g',int2str(d)],'G','R','N','d','K0','KF','Kstep','ORDmax','ORD','ResLayer')

end
