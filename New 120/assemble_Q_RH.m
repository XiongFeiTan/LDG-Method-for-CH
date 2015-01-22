function [FQ,RL] = assemble_Q_RH(U,R,M,Q_BR,Q_BL,D,N)
%function assemble_Q_RH  calculate FQ,RL
%   input£ºEnv£¬N£¬h U£¬D
%                   Q_BR,Q_BL boudary matrix 
%   output£ºFQ vector£¬RL 

%predistribution
FQ=zeros(3*N,1,'double');
RL=zeros(3*N,1,'double');

%boundary processing
for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
end
RL(1:3,1)=R(3*N-2:3*N,1);

%right hand item
FQ=M*U+(D-Q_BR)*R+Q_BL*RL;
end