function mainDriver
clear;clc;
%%
N=640;a=0;b=30;
%degree of freedom
Dof=3;
%The total degree of freedom
tDof=Dof*N;
%% 
T=1000;
dt=1.0e-2;
NT=round(T/dt);
%% generate mesh coordinate value
cv=zeros(1,(N+1));
cv=linspace(a,b,N+1);
%the gauss points and its weights
p=[-0.9602898565 -0.7966664774 -0.5255324099 -0.1834346425 ...
    0.1834346425 0.5255324099 0.7966664774 0.9602898565]';
w=[0.1012285363 0.2223810345 0.3137066459 0.3626837834 ...
   0.3626837834 0.3137066459 0.2223810345 0.1012285363]';
%% compute the integral  on a element & mapping value & jacobi weight 
[mv, jw]=elemInt(cv(1), cv(2), p, w);
% calculate the element mass matrix
[phi1,phi2,phi3]=basisfun(cv(1), cv(2), mv,'Mbasis');
m=zeros(Dof,Dof);
m(1,1)=sum(phi1.*phi1.*jw);
m(1,2)=sum(phi2.*phi1.*jw);
m(1,3)=sum(phi3.*phi1.*jw);
%
m(2,1)=sum(phi1.*phi2.*jw);
m(2,2)=sum(phi2.*phi2.*jw);
m(2,3)=sum(phi3.*phi2.*jw);
%
m(3,1)=sum(phi1.*phi3.*jw);
m(3,2)=sum(phi2.*phi3.*jw);
m(3,3)=sum(phi3.*phi3.*jw);
% calculate the element d mass matrix
[dphi1,dphi2,dphi3]=basisfun(cv(1), cv(2), mv,'Dbasis');
% calculate the element d mass matrix
 d=zeros(Dof,Dof);
d(1,1)=sum(phi1.*dphi1.*jw);
d(1,2)=sum(phi2.*dphi1.*jw);
d(1,3)=sum(phi3.*dphi1.*jw);
%
d(2,1)=sum(phi1.*dphi2.*jw);
d(2,2)=sum(phi2.*dphi2.*jw);
d(2,3)=sum(phi3.*dphi2.*jw);
 %
d(3,1)=sum(phi1.*dphi3.*jw);
d(3,2)=sum(phi2.*dphi3.*jw);
d(3,3)=sum(phi3.*dphi3.*jw);
%%  calculate U0 and predeclare the space for storage the numerical solution
U=zeros(tDof,1);
% discrete the initial value by piecewise L2 projection,right hand side
for j=1:N
    [rhs] = initialzieU(Dof,cv(j),cv(j+1),p,w);
	U(Dof*(j-1)+1:Dof*(j-1)+Dof)=m\rhs;
end
% U0 done
%% calculate R and predeclare the space for storage the numerical solution
 R=zeros(tDof,1);
for j=1:N
    uphi=U(3*(j-1)+1)*phi1+U(3*(j-1)+2)*phi2+U(3*(j-1)+3)*phi3;
    rhs=zeros(3,1);
    rhs(1)=-sum(dphi1.*uphi.*jw);
    rhs(2)=-sum(dphi2.*uphi.*jw);
    rhs(3)=-sum(dphi3.*uphi.*jw); 
    if(j==1)
        rrur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        rrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        rlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        rlul=U(3*(N-1)+1)+U(3*(N-1)+2)+U(3*(N-1)+3);
        %
        rhs(1)=rhs(1)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
        rhs(2)=rhs(2)+(rrur+rrul)/2.0+(rlur+rlul)/2.0;
        rhs(3)=rhs(3)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
    elseif(j==N)
        rrur=U(1)-U(2)+U(3);
        rrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        rlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        rlul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
        %
        rhs(1)=rhs(1)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
        rhs(2)=rhs(2)+(rrur+rrul)/2.0+(rlur+rlul)/2.0;
        rhs(3)=rhs(3)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
    else
        rrur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        rrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        rlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        rlul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
       %
        rhs(1)=rhs(1)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
        rhs(2)=rhs(2)+(rrur+rrul)/2.0+(rlur+rlul)/2.0;
        rhs(3)=rhs(3)+(rrur+rrul)/2.0-(rlur+rlul)/2.0;
    end
	R(Dof*(j-1)+1:Dof*(j-1)+Dof)=m\rhs;
end
% R done
%% calculate Q and predeclare the space for storage the numerical solution
Q=zeros(tDof,1);
for j=1:N
    uphi=U(3*(j-1)+1)*phi1+U(3*(j-1)+2)*phi2+U(3*(j-1)+3)*phi3;
    rphi=R(3*(j-1)+1)*phi1+R(3*(j-1)+2)*phi2+R(3*(j-1)+3)*phi3;
    rhs=zeros(3,1);
    rhs(1)=sum(dphi1.*rphi.*jw)+sum(phi1.*uphi.*jw);
    rhs(2)=sum(dphi2.*rphi.*jw)+sum(phi2.*uphi.*jw);
    rhs(3)=sum(dphi3.*rphi.*jw)+sum(phi3.*uphi.*jw); 
    if(j==1)
        qrrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        qrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrl=R(3*(N-1)+1)+R(3*(N-1)+2)+R(3*(N-1)+3);
        %
        rhs(1)=rhs(1)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
        rhs(2)=rhs(2)-(qrrr+qrrl)/2.0-(qlrr+qlrl)/2.0;
        rhs(3)=rhs(3)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
    elseif(j==N)
        qrrr=R(1)-R(2)+R(3);
        qrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        rhs(1)=rhs(1)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
        rhs(2)=rhs(2)-(qrrr+qrrl)/2.0-(qlrr+qlrl)/2.0;
        rhs(3)=rhs(3)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
    else
        qrrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        qrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qlrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        rhs(1)=rhs(1)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
        rhs(2)=rhs(2)-(qrrr+qrrl)/2.0-(qlrr+qlrl)/2.0;
        rhs(3)=rhs(3)-(qrrr+qrrl)/2.0+(qlrr+qlrl)/2.0;
    end
	Q(Dof*(j-1)+1:Dof*(j-1)+Dof)=m\rhs;
end
% Q done 
%% calculate P and predeclare the space for storage the numerical solution
P=zeros(tDof,1);
for j=1:N
    uphi=U(3*(j-1)+1)*phi1+U(3*(j-1)+2)*phi2+U(3*(j-1)+3)*phi3;
    rphi=R(3*(j-1)+1)*phi1+R(3*(j-1)+2)*phi2+R(3*(j-1)+3)*phi3;
    rhs=zeros(3,1);
    urphi=uphi.*rphi;
    rhs(1)=-sum(dphi1.*urphi.*jw);
    rhs(2)=-sum(dphi2.*urphi.*jw);
    rhs(3)=-sum(dphi3.*urphi.*jw); 
    if(j==1)
        prrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        prrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        plrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        plrl=R(3*(N-1)+1)+R(3*(N-1)+2)+R(3*(N-1)+3);
        %
        prur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        prul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        plur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        plul=U(3*(N-1)+1)+U(3*(N-1)+2)+U(3*(N-1)+3);
        %
        rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
    elseif(j==N)
        prrr=R(1)-R(2)+R(3);
        prrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        plrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        plrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        prur=U(1)-U(2)+U(3);
        prul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        plur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        plul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
        %
        rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;      
    else
        prrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        prrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        plrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        plrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        prur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        prul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        plur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        plul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
        %
        rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
    end
	P(Dof*(j-1)+1:Dof*(j-1)+Dof)=m\rhs;
end
% P done
%%  calculate QT and predeclare the space for storage the numerical solution
QT=zeros(tDof,1);
for j=1:N  
    uphi=U(3*(j-1)+1)*phi1+U(3*(j-1)+2)*phi2+U(3*(j-1)+3)*phi3;
    rphi=R(3*(j-1)+1)*phi1+R(3*(j-1)+2)*phi2+R(3*(j-1)+3)*phi3;
    pphi=P(3*(j-1)+1)*phi1+P(3*(j-1)+2)*phi2+P(3*(j-1)+3)*phi3;
    f=3/2.0.*(uphi.*uphi)-pphi+1/2.0.*(rphi.*rphi);
    rhs=zeros(3,1);
    rhs(1)=sum(dphi1.*f.*jw);
    rhs(2)=sum(dphi2.*f.*jw);
    rhs(3)=sum(dphi3.*f.*jw); 
    if(j==1)
        qtrrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        qtrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrl=R(3*(N-1)+1)+R(3*(N-1)+2)+R(3*(N-1)+3); 
        %
        qtrur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        qtrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlul=U(3*(N-1)+1)+U(3*(N-1)+2)+U(3*(N-1)+3);
        %
        qtrpr=P(3*j+1)-P(3*j+2)+P(3*j+3);
        qtrpl=P(3*(j-1)+1)+P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpr=P(3*(j-1)+1)-P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpl=P(3*(N-1)+1)+P(3*(N-1)+2)+P(3*(N-1)+3);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
    elseif(j==N)
        qtrrr=R(1)-R(2)+R(3);
        qtrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        qtrur=U(1)-U(2)+U(3);
        qtrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
        %
        qtrpr=P(1)-P(2)+P(3);
        qtrpl=P(3*(j-1)+1)+P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpr=P(3*(j-1)+1)-P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpl=P(3*(j-2)+1)+P(3*(j-2)+2)+P(3*(j-2)+3);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
    else
        qtrrr=R(3*j+1)-R(3*j+2)+R(3*j+3);
        qtrrl=R(3*(j-1)+1)+R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrr=R(3*(j-1)+1)-R(3*(j-1)+2)+R(3*(j-1)+3);
        qtlrl=R(3*(j-2)+1)+R(3*(j-2)+2)+R(3*(j-2)+3);
        %
        qtrur=U(3*j+1)-U(3*j+2)+U(3*j+3);
        qtrul=U(3*(j-1)+1)+U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlur=U(3*(j-1)+1)-U(3*(j-1)+2)+U(3*(j-1)+3);
        qtlul=U(3*(j-2)+1)+U(3*(j-2)+2)+U(3*(j-2)+3);
        %
        qtrpr=P(3*j+1)-P(3*j+2)+P(3*j+3);
        qtrpl=P(3*(j-1)+1)+P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpr=P(3*(j-1)+1)-P(3*(j-1)+2)+P(3*(j-1)+3);
        qtlpl=P(3*(j-2)+1)+P(3*(j-2)+2)+P(3*(j-2)+3);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
    end
    QT(Dof*(j-1)+1:Dof*(j-1)+Dof)=m\rhs;
end
% QT calcuate done
%% coef of coupling matrix
jr=[-1/2,1/2,-1/2;-1/2,1/2,-1/2;-1/2,1/2,-1/2];
j=[0,-1,0;-1,0,-1;0,-1,0];
jl=[1/2,1/2,1/2;-1/2,-1/2,-1/2;1/2,1/2,1/2];
%%  coupling matrix & RK method
[C,M]= assemble(Dof,N,m,d,jr,j,jl);
[Q] = RK3(tDof,Q,QT,dt,M,C,m,N,jw,phi1,phi2,phi3,dphi1,dphi2,dphi3);
%%  time-discrete
for K=1:NT
    [QT1]  = getQT(tDof,Q,M,C,m,N,jw,phi1,phi2,phi3,dphi1,dphi2,dphi3);
    [Q] = RK3(tDof,Q,QT1,dt,M,C,m,N,jw,phi1,phi2,phi3,dphi1,dphi2,dphi3);
end
%%  calculate U
U=zeros(tDof,1);
R=zeros(tDof,1);
F=zeros(tDof*2,1);
URV=zeros(tDof*2,1);
%%
F(1:tDof,1)=M*Q;
URV=C\F;
U=URV(1:tDof,1);
R=URV(tDof+1:tDof*2,1);
%% numerical solution nu  
nu=zeros(N,1);
nu=U(1:3:size(U,1))-U(3:3:size(U,1))/2.0;
%% example3 plot one peakon exact sloution and numerical sloution
eu=zeros(N,1);
for i= 1:N
    cv(i)=(cv(i)+cv(i+1))/2;
%     eu(i,1)=exp(-abs(cv(i)-T-10));
    eu(i,1)=exp(-abs(cv(i)-T+980));
end
%--------------------exact&numerical--------------
c=linspace(a,b,N);
plot(c,eu, '-k+', c, nu, '-ro');
legend ( 'Exact','Numerical',0) ;
axis([a b -0.2,1.2]);
title ( sprintf ( 'N= %d,T=%d\n ',N, T ) );
xlabel ( 'X' );
ylabel ( 'U' );
%% example 4





%% example 5



%% example 6



%% example 7


end
