function [FR,UR] = assemble_R_RH(U,R_BR,R_BL,D,N)
%function assemble_R_RH assemble R equation's matrix
%   input£ºEnv£¬N£¬h U£¬D
%                   R_BR,R_BL boudary matrix 
%   output£ºM,D £¬FR vector£¬UR U

%predistribution
FR=zeros(3*N,1,'double');
UR=zeros(3*N,1,'double');
%boudary item processing
UR(3*N-2:3*N,1)=U(1:3,1);
for i=1:N-1
    UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
end
%right item
FR=-(R_BL+D)*U+R_BR*UR;
end
