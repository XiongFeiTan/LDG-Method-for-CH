function [Q,RR] = assemble_Q_RH(U,R,M,Q_BR,Q_BL,D,N)
%函数 assemble_Q_Elem  求Q向量
%   输入变量：N为单元数,h为步长,D为求导后的基函数做内积的总刚矩阵,M为基函数做内积的总刚矩阵
%                    R向量，U向量， Q_BR,Q_BL关于Q方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：Q向量，RR为由R变化而来，边值问题处理右端
%                   

%预分配
F=zeros(3*N,1);
RR=zeros(3*N,1);

% 边界处理
for i=1:N-1
    RR(3*i-2:3*i,1)=R(3*i+1:3*i+3,1);
end
RR(3*N-2:3*N,1)=R(1:3,1);

%右端项
F=M*U+(D+Q_BL)*R-(Q_BR*RR);
%[l,u]=lu(M);
%求线性方程
Q=M\F;
end