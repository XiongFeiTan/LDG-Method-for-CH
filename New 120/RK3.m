function [QT]=RK3(M,FQT,dt,Q,C,N,Env,h)
%function RK3 time discretization
%inputs£º
%outputs:
%%
%predistribution
%%
[l,u]=lu(M);
rk1=u\l\FQT;
%%
Q1=Q+0.5*dt*rk1;
CRH=zeros(6*N,1);
CRH(1:3*N,1)=M*Q1;
%calculate U,R
X=C\CRH;
U=X(1:3*N,1);
R=X(3*N+1:6*N,1);
%new UR,RLneed new processing U,R translation of a unit vector
UR=zeros(3*N,1);
UR(3*N-2:3*N,1)=U(1:3,1);
for i=1:N-1
        UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
end   
RL=zeros(3*N,1);
for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
end
RL(1:3,1)=R(3*N-2:3*N,1);
%re-sloved P
[P_BR,P_BL,RU,FP1] = assemble_P(U,R,UR,RL,Env,N,h);
P=M\FP1;
%re-sloved QT
[QT_BR,QT_BL,RUP,FQT1] = assemble_QT(U,R,UR,RL,P,Env,N,h);
rk2=u\l\FQT1;
%%
Q2=Q-dt*rk1+2*dt*rk2;
CRH=zeros(6*N,1);
CRH(1:3*N,1)=M*Q2;
%calculate U,R
X=C\CRH;
U=X(1:3*N,1);
R=X(3*N+1:6*N,1);
%new UR,RLneed new processing U,R translation of a unit vector
UR=zeros(3*N,1);
UR(3*N-2:3*N,1)=U(1:3,1);
for i=1:N-1
        UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
end 
RL=zeros(3*N,1);
for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
end
RL(1:3,1)=R(3*N-2:3*N,1);
 %re-sloved P
[P_BR,P_BL,RU,FP2] = assemble_P(U,R,UR,RL,Env,N,h);
P=M\FP2;
%re-sloved QT
[QT_BR,QT_BL,RUP,FQT2] = assemble_QT(U,R,UR,RL,P,Env,N,h);
rk3=u\l\FQT2;
%%
QT=Q+(dt/6)*(rk1+4*rk2+rk3);
end