function [QT] = assemble_QT_RH(M,Q,QT_BR,QT_BL,RUP,N,dt)
%函数 assemble_QT_Elem 求QT
%   输入变量：Q_BR,Q_BL关于Q方程边界右端点的总刚矩阵,左端点的总刚矩阵，RUP
%   输出变量：QT向量
%   欧拉方法                 

%预分配
F=zeros(3*N,1);
%右端项
F=M*Q+dt*(RUP+QT_BL-QT_BR);
%[l,u]=lu(M);
%求解
QT=M\F;
end