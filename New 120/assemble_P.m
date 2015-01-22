function [P_BR,P_BL,RU,FP] = assemble_P(U,R,UR,RL,Env,N,h)
%function assemble_P assemble P equation's matrix
%   input£ºEnv£¬N£¬h,UR,RL
%   output£ºRU FP  P_BR,P_BL,...boudary matrix                  
%predistribution

bl = zeros(3,1);
br = zeros(3,1);
P_BR=zeros(3*N,1,'double');
P_BL=zeros(3*N,1,'double');
RU=zeros(3*N,1,'double');
FP=zeros(3*N,1,'double');

%calculate RU£¬and assemble
for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    %three coefficient
    r(1:3,1)=R(3*i-2:3*i,1);
    u(1:3,1)=U(3*i-2:3*i,1);
    
    %integrand R£¬U£¬PHI
    F2=@(x)((r(1,1)+r(2,1)*(x-mid)+r(3,1)*(x-mid).^2).*(u(1,1)+u(2,1)*(x-mid)+u(3,1)*(x-mid).^2));
    F3=@(x)((r(1,1)+r(2,1)*(x-mid)+r(3,1)*(x-mid).^2).*(u(1,1)+u(2,1)*(x-mid)+u(3,1)*(x-mid).^2).*(2*(x-mid)));
    %quadrature
    quad1 = 0;
    quad2 = quadrature(F2,a,b);
    quad3 = quadrature(F3,a,b);
    % assemble RU
    RU(3*i-2:3*i,1)=[quad1,quad2,quad3]';
end

%boudary processing
RR=zeros(3*N,1,'double');
RR(3*N-2:3*N,1)=R(1:3,1);
for i=1:N-1
    RR(3*i-2:3*i,1)=R(3*i+1:3*i+3,1);
end

% calculate P_BR,P_BL  matrix  and assemble 
for i=1:N
    %three coefficient 
    r(1:3,1)=R(3*i-2:3*i,1);
    u(1:3,1)=U(3*i-2:3*i,1);
    ur(1:3,1)=UR(3*i-2:3*i,1);
    rl(1:3,1)=RL(3*i-2:3*i,1);
    rr(1:3,1)=RR(3*i-2:3*i,1);
    
    brr=(1/2*(rr(1,1)+rr(2,1)*(-h/2)+rr(3,1)*(h^2/4))^2-1/2*(r(1,1)+r(2,1)*(h/2)+r(3,1)*(h^2/4))^2)/((rr(1,1)+rr(2,1)*(-h/2)+rr(3,1)*(h^2/4))-(r(1,1)+r(2,1)*(h/2)+r(3,1)*(h^2/4)));
    brl=(1/2*(r(1,1)+r(2,1)*(-h/2)+r(3,1)*(h^2/4))^2-1/2*(rl(1,1)+rl(2,1)*(h/2)+rl(3,1)*(h^2/4))^2)/((r(1,1)+r(2,1)*(-h/2)+r(3,1)*(h^2/4))-(rl(1,1)+rl(2,1)*(h/2)+rl(3,1)*(h^2/4)));
    br=[brr*(ur(1,1)+ur(2,1)*(-h/2)+ur(3,1)*(h^2/4)),...
          brr*(ur(1,1)+ur(2,1)*(-h/2)+ur(3,1)*(h^2/4))*(h/2),...
          brr*(ur(1,1)+ur(2,1)*(-h/2)+ur(3,1)*(h^2/4))*(h^2/4)]';
    bl=[brl*(u(1,1)+u(2,1)*(-h/2)+u(3,1)*(h^2/4)),...
          brl*(u(1,1)+u(2,1)*(-h/2)+u(3,1)*(h^2/4))*(-h/2),...
          brl*(u(1,1)+u(2,1)*(-h/2)+u(3,1)*(h^2/4))*(h^2/4)]';

%assemble 
    P_BR(3*i-2:3*i,1)=br;
    P_BL(3*i-2:3*i,1)=bl;
end
%right hand item 
FP=-RU+P_BR-P_BL;
end


