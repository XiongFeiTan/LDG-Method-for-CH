function [QT_BR,QT_BL,RUP,FQT] = assemble_QT(U,R,UR,RL,P,Env,N,h)
%function assemble_QT QT equation's matrix
%   input£ºEnv£¬N£¬h,U,R,UR,RL,P
%   
%   output£ºRUP is u,r,p multiplication vector
%                    QT_BR,QT_BL boudary matrix
%predistribution
QT_BR=zeros(3*N,1,'double');
QT_BL=zeros(3*N,1,'double');
RUP=zeros(3*N,1,'double');
FQT=zeros(3*N,1,'double');

%calculate f,flux phi 
umax=zeros(1,N,'double');
for i=1:N
    un(1:3,1)=U(3*i-2:3*i,1);
    u1=un(1,1)+un(2,1)*(-h/2)+un(3,1)*(h^2/4);
    u2=un(1,1)+un(2,1)*(h/2)+un(3,1)*(h^2/4);
    umax(1,i)=max(u1,u2);
end
phi=max(abs(3.*umax));

%calculate RUP
for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    %three coefficient
    r(1:3,1)=R(3*i-2:3*i,1);
    u(1:3,1)=U(3*i-2:3*i,1);
    p(1:3,1)=P(3*i-2:3*i,1);
    %integrand R£¬U£¬PHI
     F2=@(x)(1/2*(r(1,1)+r(2,1)*(x-mid)+r(3,1)*(x-mid).^2).^2+3/2*(u(1,1)+u(2,1)*(x-mid)+u(3,1)*(x-mid).^2).^2-(p(1,1)+p(2,1)*(x-mid)+p(3,1)*(x-mid).^2));
    F3=@(x)((1/2*(r(1,1)+r(2,1)*(x-mid)+r(3,1)*(x-mid).^2).^2+3/2*(u(1,1)+u(2,1)*(x-mid)+u(3,1)*(x-mid).^2).^2-(p(1,1)+p(2,1)*(x-mid)+p(3,1)*(x-mid).^2)).*(2*(x-mid)));
    %quadrature
    quad1 = 0;
    quad2 = quadrature(F2,a,b);
    quad3 = quadrature(F3,a,b);
    % assemble RUP
    RUP(3*i-2:3*i,1)=[quad1,quad2,quad3]';
end
    
%boudary item processing
PL=zeros(3*N,1,'double');
for i=2:N
    PL(3*i-2:3*i,1)=P(3*i-5:3*i-3,1);
end
PL(1:3,1)=P(3*N-2:3*N,1);

%boudary item processing
UL=zeros(3*N,1,'double');
for i=2:N
    UL(3*i-2:3*i,1)=U(3*i-5:3*i-3,1);
end
UL(1:3,1)=U(3*N-2:3*N,1);


%element loop assemble
for i=1:N
    %three coefficient
    r(1:3,1)=R(3*i-2:3*i,1);
    u(1:3,1)=U(3*i-2:3*i,1);
    p(1:3,1)=P(3*i-2:3*i,1);
    pl(1:3,1)=PL(3*i-2:3*i,1);
    ur(1:3,1)=UR(3*i-2:3*i,1);
    ul(1:3,1)=UL(3*i-2:3*i,1);
    rl(1:3,1)=RL(3*i-2:3*i,1);
  %f flux value
  fur=1/2*(3/2*(ur(1,1)+ur(2,1)*(-h/2)+ur(3,1)*(h^2/4))^2+3/2*(u(1,1)+u(2,1)*(h/2)+u(3,1)*(h^2/4))^2-phi*((ur(1,1)+ur(2,1)*(-h/2)+ur(3,1)*(h^2/4))-(u(1,1)+u(2,1)*(h/2)+u(3,1)*(h^2/4))));
  ful=1/2*(3/2*(u(1,1)+u(2,1)*(-h/2)+u(3,1)*(h^2/4))^2+3/2*(ul(1,1)+ul(2,1)*(h/2)+ul(3,1)*(h^2/4))^2-phi*((u(1,1)+u(2,1)*(-h/2)+u(3,1)*(h^2/4))-(ul(1,1)+ul(2,1)*(h/2)+ul(3,1)*(h^2/4))));
    
%
    br1=(1/2*(r(1,1)+r(2,1)*(h/2)+r(3,1)*(h^2/4))^2+fur-(p(1,1)+p(2,1)*(h/2)+p(3,1)*(h^2/4)));
    br2=(1/2*(r(1,1)+r(2,1)*(h/2)+r(3,1)*(h^2/4))^2+fur-(p(1,1)+p(2,1)*(h/2)+p(3,1)*(h^2/4)))*(h/2);
    br3=(1/2*(r(1,1)+r(2,1)*(h/2)+r(3,1)*(h^2/4))^2+fur-(p(1,1)+p(2,1)*(h/2)+p(3,1)*(h^2/4)))*(h^2/4);
    br=[br1,br2,br3]';
    %assemble
   QT_BR(3*i-2:3*i,1)=br;
    
    bl1=(1/2*(rl(1,1)+rl(2,1)*(h/2)+rl(3,1)*(h^2/4))^2+ful-(pl(1,1)+pl(2,1)*(h/2)+pl(3,1)*(h^2/4)));
    bl2=(1/2*(rl(1,1)+rl(2,1)*(h/2)+rl(3,1)*(h^2/4))^2+ful-(pl(1,1)+pl(2,1)*(h/2)+pl(3,1)*(h^2/4)))*(-h/2);
    bl3=(1/2*(rl(1,1)+rl(2,1)*(h/2)+rl(3,1)*(h^2/4))^2+ful-(pl(1,1)+pl(2,1)*(h/2)+pl(3,1)*(h^2/4)))*(h^2/4);
    
     bl=[bl1,bl2,bl3]';
     %assemble 
    QT_BL(3*i-2:3*i,1)=bl;
end
%right hand item
FQT=RUP+QT_BL-QT_BR;
end


