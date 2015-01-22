function [U] = Initializeu(Env,M,N)
%function  Initializeu calcualte the initial value U
%   input£ºEnv£¬N element number £¬M inner product of the total mass basis function matrix
%   output£ºU

%predistribution
U=zeros(3*N,1,'double');
%right hand item 
F=zeros(3*N,1,'double');
% element loop 
for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    % integrand
    f1=@(x)(0.25*exp(-abs(x)));
    f2=@(x)(0.25*exp(-abs(x)).*(x-mid));
    f3=@(x)(0.25*exp(-abs(x)).*(x-mid).^2);
    % quadrature
    quad1=quadrature(f1,a,b);
    quad2=quadrature(f2,a,b);
    quad3=quadrature(f3,a,b);
    %assemble right item
    F(3*i-2:3*i)=[quad1,quad2,quad3]';
end
%[l,u]=lu(M);
%slove the system of linear equations 
U=M\F;
end


