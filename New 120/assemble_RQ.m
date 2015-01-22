function [r_br,r_bl,q_br,q_bl,M,D,R_BR,R_BL,Q_BR,Q_BL] = assemble_RQ(Env,N,h)
%function assemble_RQ assemble R,Q equation's matrix
%   input£ºEnv£¬N£¬h
%   output£ºM,D
%                    R_BR,R_BL,...boudary matrix 
%                    r_br,r_bl,...local boudary matrix
%predistribution
m = zeros(3,3);
d = zeros(3,3);
r_bl = zeros(3,3);
r_br = zeros(3,3);
q_bl = zeros(3,3);
q_br = zeros(3,3);
% matrix
M=zeros(3*N,3*N,'double');
D=zeros(3*N,3*N,'double');
R_BR=zeros(3*N,3*N,'double');
R_BL=zeros(3*N,3*N,'double');
Q_BR=zeros(3*N,3*N,'double');
Q_BL=zeros(3*N,3*N,'double');

%local matrix
r_br=[1,-h/2,h^2/4;h/2,-h^2/4,h^3/8;h^2/4,-h^3/8,h^4/16];
r_bl=[1,-h/2,h^2/4;-h/2,h^2/4,-h^3/8;h^2/4,-h^3/8,h^4/16];
q_br=[1,h/2,h^2/4;h/2,h^2/4,h^3/8;h^2/4,h^3/8,h^4/16];
q_bl=[1,h/2,h^2/4;-h/2,-h^2/4,-h^3/8;h^2/4,h^3/8,h^4/16];

%basis function
% phi1=@(x)(1);
% phi2=@(x)(x-mid);
% phi3=@(x)((x-mid)^2);
% derivation basis function
% Dphi1=@(x)(0);
% Dphi2=@(x)(1);
% Dphi3=@(x)(2*(x-mid));
%element loop 
for i=1:N
    %a left endpoint 
    a=Env(i,1);
    %b right endpoint
    b=Env(i,2);
    mid=(a+b)/2;
    %basis function multiplication
    f1=@(x)(1);
    f2=@(x)(x-mid);
    f3=@(x)((x-mid).^2);
    f4=@(x)((x-mid).^3);
    f5=@(x)((x-mid).^4);

    %derivation basis function multiplication 
    y1=@(x)(1);
    y2=@(x)(x-mid);
    y3=@(x)((x-mid).^2);
    y4=@(x)(2*(x-mid));
    y5=@(x)(2*(x-mid).^2);
    y6=@(x)(2*(x-mid).^3);
    %quadrature 
    quad1 = quadrature(f1,a,b);
    quad2 = quadrature(f2,a,b);
    quad3 = quadrature(f3,a,b);
    quad4 = quadrature(f4,a,b);
    quad5 = quadrature(f5,a,b);
    %
    dquad1 = quadrature(y1,a,b);
    dquad2 = quadrature(y2,a,b);
    dquad3 = quadrature(y3,a,b);
    dquad4 = quadrature(y4,a,b);
    dquad5 = quadrature(y5,a,b);
    dquad6 = quadrature(y6,a,b);
    % m matrix 
    m=[quad1,quad2,quad3;quad2,quad3,quad4;quad3,quad4,quad5];
    % d matrix
    d=[0,0,0;dquad1,dquad2,dquad3;dquad4,dquad5,dquad6];
    % assemble
    M(3*i-2:3*i,3*i-2:3*i)=m;
    D(3*i-2:3*i,3*i-2:3*i)=d;
    R_BR(3*i-2:3*i,3*i-2:3*i)=r_br;
    R_BL(3*i-2:3*i,3*i-2:3*i)=r_bl;
    Q_BR(3*i-2:3*i,3*i-2:3*i)=q_br;
    Q_BL(3*i-2:3*i,3*i-2:3*i)=q_bl;
    end
end


