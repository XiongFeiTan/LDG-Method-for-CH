function [FQ,RL] = assemble_Q_RH(U,R,M,Q_BR,Q_BL,D,N)
%函数 assemble_Q_RH 求FQ向量,RL
%   输入变量：N为单元数,h为步长,D为求导后的基函数做内积的总刚矩阵,M为基函数做内积的总刚矩阵
%                    R向量，U向量， Q_BR,Q_BL关于Q方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：FQ向量，RL为由R变化而来，边值问题处理
               

%预分配
FQ=zeros(3*N,1);
RL=zeros(3*N,1);

% 边界处理
for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
end
RL(1:3,1)=R(3*N-2:3*N,1);

%右端项
FQ=M*U+(D-Q_BR)*R+Q_BL*RL;
end