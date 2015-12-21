function mainDriver
clear;
clc;
%%
N=120;a=-25;b=25;
%degree of freedom
Dof=4;
%The total degree of freedom
tDof=Dof*N;
%% time
T=1;
dt=0.01;
NT=round(T/dt);
%%  generate mesh coordinate value
% cv=zeros(1,N+1);
cv=linspace(a,b,N+1);
%%  the gauss points and its weights
p=[-0.9602898565 -0.7966664774 -0.5255324099 -0.1834346425 ...
    0.1834346425 0.5255324099 0.7966664774 0.9602898565]';
w=[0.1012285363 0.2223810345 0.3137066459 0.3626837834 ...
   0.3626837834 0.3137066459 0.2223810345 0.1012285363]';
%%  compute the integral  on a element % mapping value % jacobi weight
[mv, jw]=elemInt(cv(1), cv(2), p, w);
% calculate the element mass matrix
[phi1,phi2,phi3,phi4]=basisfun(cv(1), cv(2), mv,'Mbasis');
m=zeros(4,4);
m(1,1)=sum(phi1.*phi1.*jw);
m(1,2)=sum(phi2.*phi1.*jw);
m(1,3)=sum(phi3.*phi1.*jw);
m(1,4)=sum(phi4.*phi1.*jw);
%
m(2,1)=sum(phi1.*phi2.*jw);
m(2,2)=sum(phi2.*phi2.*jw);
m(2,3)=sum(phi3.*phi2.*jw);
m(2,4)=sum(phi4.*phi2.*jw);
%
m(3,1)=sum(phi1.*phi3.*jw);
m(3,2)=sum(phi2.*phi3.*jw);
m(3,3)=sum(phi3.*phi3.*jw);
m(3,4)=sum(phi4.*phi3.*jw);
%
m(4,1)=sum(phi1.*phi4.*jw);
m(4,2)=sum(phi2.*phi4.*jw);
m(4,3)=sum(phi3.*phi4.*jw);
m(4,4)=sum(phi4.*phi4.*jw);
%% calculate the derivate element mass matrix
[dphi1,dphi2,dphi3,dphi4]=basisfun(cv(1), cv(2), mv,'Dbasis');
d=zeros(4,4);
d(1,1)=sum(phi1.*dphi1.*jw);
d(1,2)=sum(phi2.*dphi1.*jw);
d(1,3)=sum(phi3.*dphi1.*jw);
d(1,4)=sum(phi4.*dphi1.*jw);
%
d(2,1)=sum(phi1.*dphi2.*jw);
d(2,2)=sum(phi2.*dphi2.*jw);
d(2,3)=sum(phi3.*dphi2.*jw);
d(2,4)=sum(phi4.*dphi2.*jw);
 %
d(3,1)=sum(phi1.*dphi3.*jw);
d(3,2)=sum(phi2.*dphi3.*jw);
d(3,3)=sum(phi3.*dphi3.*jw);
d(3,4)=sum(phi4.*dphi3.*jw);
%
d(4,1)=sum(phi1.*dphi4.*jw);
d(4,2)=sum(phi2.*dphi4.*jw);
d(4,3)=sum(phi3.*dphi4.*jw);
d(4,4)=sum(phi4.*dphi4.*jw);
%% calculate U0 & predeclare the space 
U=zeros(tDof,1);
for j=1:N
    [rhs] = initialzieU(Dof,cv(j),cv(j+1),p,w);
	U(4*(j-1)+1:4*(j-1)+4)=m\rhs;
end
% U0 done
%%  calculate R & predeclare the space
UR=[U(5:4*N);U(1:4)];
UL=[U(4*N-3:4*N);U(1:4*N-4)];
M=kron(speye(N),m);
D=kron(speye(N),d);
jl=[1/2,1/2,1/2,1/2;-1/2,-1/2,-1/2,-1/2;1/2,1/2,1/2,1/2;-1/2,-1/2,-1/2,-1/2];
j=[0,-1,0,-1;-1,0,-1,0;0,-1,0,-1;-1,0,-1,0];
jr=[-1/2,1/2,-1/2,1/2;-1/2,1/2,-1/2,1/2;-1/2,1/2,-1/2,1/2;-1/2,1/2,-1/2,1/2];
BR=kron(speye(N),jr);
B=kron(speye(N),j);
BL=kron(speye(N),jl);
R=M\(-(D+B)*U-BR*UR-BL*UL);
%% calculate Q & predeclare the space
RL=[R(4*N-3:4*N);R(1:4*N-4)];
RR=[R(5:4*N);R(1:4)];
Q=M\(M*U+(D+B)*R+BR*RR+BL*RL);
%%
P=zeros(tDof,1);
for j=1:N
    uphi=U(4*(j-1)+1)*phi1+U(4*(j-1)+2)*phi2+U(4*(j-1)+3)*phi3+U(4*(j-1)+4)*phi4;
    rphi=R(4*(j-1)+1)*phi1+R(4*(j-1)+2)*phi2+R(4*(j-1)+3)*phi3+R(4*(j-1)+4)*phi4;
    rhs=zeros(4,1);
    rhs(1)=-sum(dphi1.*rphi.*uphi.*jw);
    rhs(2)=-sum(dphi2.*rphi.*uphi.*jw);
    rhs(3)=-sum(dphi3.*rphi.*uphi.*jw); 
    rhs(4)=-sum(dphi4.*rphi.*uphi.*jw); 
    if(j==1)
        plrl=R(4*(N-1)+1)+R(4*(N-1)+2)+R(4*(N-1)+3)+R(4*(N-1)+4);
        plrr=R(4*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        prrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        prrr=R(4*j+1)-R(4*j+2)+R(4*j+3)-R(4*j+4);
        %
        plur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*(j-1)+4);
        prur=U(4*j+1)-U(4*j+2)+U(4*j+3)-U(4*j+4);
        prul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*(j-1)+4);
        plul=U(4*(N-1)+1)+U(4*(N-1)+2)+U(4*(N-1)+3)+U(4*(N-1)+4);  
        %
        rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(4)=rhs(4)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
    elseif(j==N)
        plrl=R(4*(j-2)+1)+R(4*(j-2)+2)+R(4*(j-2)+3)+R(4*(j-1)+4);
        plrr=R(4*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        prrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        prrr=R(1)-R(2)+R(3)-R(4);
        %
        plur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*(j-1)+4);
        prur=U(1)-U(2)+U(3)-U(4);
        prul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*(j-1)+4);
        plul=U(4*(j-2)+1)+U(4*(j-2)+2)+U(4*(j-2)+3)+U(4*(j-2)+4);  
        %
       rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(4)=rhs(4)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
    else
        plrl=R(4*(j-2)+1)+R(4*(j-2)+2)+R(4*(j-2)+3)+R(4*(j-1)+4);
        plrr=R(4*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        prrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        prrr=R(4*j+1)-R(4*j+2)+R(4*j+3)-R(4*j+4);
        %
        plur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*(j-1)+4);
        prur=U(4*j+1)-U(4*j+2)+U(4*j+3)-U(4*j+4);
        prul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*(j-1)+4);
        plul=U(4*(j-2)+1)+U(4*(j-2)+2)+U(4*(j-2)+3)+U(4*(j-2)+4);  
        %
        rhs(1)=rhs(1)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(2)=rhs(2)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
        rhs(3)=rhs(3)+(prrr+prrl)*(prur+prul)/4.0-(plrr+plrl)*(plur+plul)/4.0;
        rhs(4)=rhs(4)+(prrr+prrl)*(prur+prul)/4.0+(plrr+plrl)*(plur+plul)/4.0;
    end
	P(4*(j-1)+1:4*(j-1)+4)=m\rhs;
end
% P done
%%  calculate QT & predeclare the space
QT=zeros(tDof,1);
for j=1:N  
    uphi=U(4*(j-1)+1)*phi1+U(4*(j-1)+2)*phi2+U(4*(j-1)+3)*phi3+U(4*(j-1)+4)*phi4;
    rphi=R(4*(j-1)+1)*phi1+R(4*(j-1)+2)*phi2+R(4*(j-1)+3)*phi3+R(4*(j-1)+4)*phi4;
    pphi=P(4*(j-1)+1)*phi1+P(4*(j-1)+2)*phi2+P(4*(j-1)+3)*phi3+P(4*(j-1)+4)*phi4;
    f=3/2.0.*(uphi.*uphi)-pphi+(rphi.*rphi)./2.0;
    rhs=zeros(4,1);
    rhs(1)=sum(dphi1.*f.*jw);
    rhs(2)=sum(dphi2.*f.*jw);
    rhs(3)=sum(dphi3.*f.*jw); 
    rhs(4)=sum(dphi4.*f.*jw);  
    if(j==1)
        qtlrl=R(4*(N-1)+1)+R(4*(N-1)+2)+R(4*(N-1)+3)+R(4*(N-1)+4);
        qtlrr=R(4*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        qtrrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        qtrrr=R(4*j+1)-R(4*j+2)+R(4*j+3)-R(4*j+4);
        %
        qtlur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*j+4);
        qtlul=U(4*(N-1)+1)+U(4*(N-1)+2)+U(4*(N-1)+3)+U(4*(N-1)+4);
        qtrur=U(4*j+1)-U(4*j+2)+U(4*j+3)-U(4*j+4);
        qtrul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*j+4);
        %
        qtrpl=P(4*(j-1)+1)+P(4*(j-1)+2)+P(4*(j-1)+3)+P(4*(j-1)+4);
        qtlpl=P(4*(N-1)+1)+P(4*(N-1)+2)+P(4*(N-1)+3)+P(4*(N-1)+4);
        qtrpr=P(4*j+1)-P(4*j+2)+P(4*j+3)-P(4*j+4);
        qtlpr=P(4*(j-1)+1)-P(4*(j-1)+2)+P(4*(j-1)+3)-P(4*(j-1)+4);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(4)=rhs(4)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
    elseif(j==N)
        qtrrr=R(1)-R(2)+R(3)-R(4);
        qtrrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        qtlrr=R(3*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        qtlrl=R(4*(j-2)+1)+R(4*(j-2)+2)+R(4*(j-2)+3)+R(4*(j-1)+4);
        %
        qtrur=U(1)-U(2)+U(3)-U(4);
        qtrul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*(j-1)+4);
        qtlur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*(j-1)+4);
        qtlul=U(4*(j-2)+1)+U(4*(j-2)+2)+U(4*(j-2)+3)+U(4*(j-2)+4);
        %
        qtrpl=P(4*(j-1)+1)+P(4*(j-1)+2)+P(4*(j-1)+3)+P(4*(j-1)+4);
        qtlpl=P(4*(j-2)+1)+P(4*(j-2)+2)+P(4*(j-2)+3)+P(4*(j-2)+4);
        qtrpr=P(1)-P(2)+P(3)-P(4);
        qtlpr=P(4*(j-1)+1)-P(4*(j-1)+2)+P(4*(j-1)+3)-P(4*(j-1)+4);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(4)=rhs(4)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
    else
        qtrrr=R(4*j+1)-R(4*j+2)+R(4*j+3)-R(4*j+4);
        qtrrl=R(4*(j-1)+1)+R(4*(j-1)+2)+R(4*(j-1)+3)+R(4*(j-1)+4);
        qtlrr=R(4*(j-1)+1)-R(4*(j-1)+2)+R(4*(j-1)+3)-R(4*(j-1)+4);
        qtlrl=R(4*(j-2)+1)+R(4*(j-2)+2)+R(4*(j-2)+3)+R(4*(j-2)+4);
        %
        qtrur=U(4*j+1)-U(4*j+2)+U(4*j+3)-U(4*j+4);
        qtrul=U(4*(j-1)+1)+U(4*(j-1)+2)+U(4*(j-1)+3)+U(4*(j-1)+4);
        qtlur=U(4*(j-1)+1)-U(4*(j-1)+2)+U(4*(j-1)+3)-U(4*(j-1)+4);
        qtlul=U(4*(j-2)+1)+U(4*(j-2)+2)+U(4*(j-2)+3)+U(4*(j-2)+4);
        %
        qtrpl=P(4*(j-1)+1)+P(4*(j-1)+2)+P(4*(j-1)+3)+P(4*(j-1)+4);
        qtlpl=P(4*(j-2)+1)+P(4*(j-2)+2)+P(4*(j-2)+3)+P(4*(j-2)+4);
        qtrpr=P(4*j+1)-P(4*j+2)+P(4*j+3)-P(4*j+4);
        qtlpr=P(4*(j-1)+1)-P(4*(j-1)+2)+P(4*(j-1)+3)-P(4*(j-1)+4);
        %
        fr=(qtrur*qtrur+qtrur*qtrul+qtrul*qtrul)/2.0;
        fl=(qtlur*qtlur+qtlur*qtlul+qtlul*qtlul)/2.0;
        Br=(qtrrr*qtrrr+qtrrl*qtrrl)/4.0;
        Bl=(qtlrr*qtlrr+qtlrl*qtlrl)/4.0;
        %
        rhs(1)=rhs(1)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(2)=rhs(2)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(3)=rhs(3)-(fr+Br-(qtrpr+qtrpl)/2.0)+(fl+Bl-(qtlpr+qtlpl)/2.0);
        rhs(4)=rhs(4)-(fr+Br-(qtrpr+qtrpl)/2.0)-(fl+Bl-(qtlpr+qtlpl)/2.0);
    end
	QT(4*(j-1)+1:4*(j-1)+4)=m\rhs;
end
% QT  done
%%  coupling matrix & RK method &  calculate Q1
[C,M]= assemble(Dof,N,m,d,jr,j,jl);
[Q1] = RK3(tDof,Q,QT,dt,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4);
%% time discrete
for K=1:NT
    [QT1]  = getQT(tDof,Q1,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4);
    [Q1] = RK3(tDof,Q1,QT1,dt,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4);
end
%% calculate the U from Q
F=sparse(tDof*2,1);
F(1:tDof,1)=M*Q1;
% URV=zeros(2*tDof,1);
URV=C\F;
U=URV(1:tDof,1);
% R=URV(tDof+1:2*tDof,1);
%% the numerical solution % exact solution
% nu=zeros(N,1);
nu=U(1:4:size(U,1))-U(3:4:size(U,1))/2.0
%% plot example2 exact sloution and numerical sloution
eu=zeros(N,1);
for i= 1:N
    cv(i)=(cv(i)+cv(i+1))/2;
    if T<=100
        eu(i,1)=0.25*exp(-abs(cv(i)-0.25*T));
    else
        eu(i,1)=0.25*exp(-abs(cv(i)-0.25*T+50));
    end
end
%------------------------exact&numerical-------------------- 
c=linspace(a,b,N);
plot(c,eu, '-k+', c, nu, '-ro');
legend ( 'Exact','Numerical',0) ;
u_min=min(min(nu));
u_max=max(max(nu));
axis([a b u_min-0.2,u_max+0.2]);
title ( sprintf ( 'N= %d,T=%d\n ',N, T ) );
xlabel ( 'X' ) ;
ylabel ( 'U' );
%%  example3 plot one peakon exact sloution and numerical sloution
% eu=zeros(N,1);
% for i= 1:N
%     cv(i)=(cv(i)+cv(i+1))/2;
%     eu(i,1)=exp(-abs(cv(i)-T+50));
% %   eu(i,1)=exp(-abs(cv(i)-10));
% end
% %--------------------exact&numerical--------------
% c=linspace(a,b,N);
% plot(c,eu, '-k+', c, nu, '-ro');
% legend ( 'Exact','Numerical',0) ;
% %u_min=min(min(nu));
% %u_max=max(max(nu));
% %axis([a b u_min-0.2,u_max+0.2]);
% axis([a b -0.2,1.2]);
% title ( sprintf ( 'N= %d,T=%d\n ',N, T ) );
% xlabel ( 'X' ) ;
% ylabel ( 'U' );
end