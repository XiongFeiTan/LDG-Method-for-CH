function [U] = Calculate_U(M,Q,D,R,Q_BL,RL,Q_BR,N)
%function Calculate Q equation's U
%   input £ºN,h,Q_BR,Q_BL about Q equation boudary matrix
%                   M,D
%  output £º    U          

%predistribution
F=zeros(3*N,1,'double');
U=zeros(3*N,1,'double');
%assemble right hand item 
F=M*Q-(D-Q_BR)*R-Q_BL*RL;
%slove 
U=M\F;
end