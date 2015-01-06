function [U] = Calculate_U(M,Q,D,R,Q_BL,RR,Q_BR,N,h)
%function Calculate_U
%   output variable
%   input variable

%predistribution
F=zeros(3*N,1);
% % impose boundary
% for i=1:N-1
%     RR(3*i-2:3*i,1)=R(3*i+1:3*i+3,1);
% end
% RR(3*N-2:3*N,1)=R(1:3,1);

%right hand side assemble
F=M*Q-(D+Q_BL)*R+Q_BR*RR;
%[l,u]=lu(M);
U=M\F;
end