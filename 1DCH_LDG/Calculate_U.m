function [U] = Calculate_U(M,Q,D,R,Q_BL,RL,Q_BR,N)
%函数Calculate Q方程里的U
%   输入变量：N为单元数，h为步长,Q_BR,Q_BL关于Q方程边界右端点的总刚矩阵,左端点的总刚矩阵
%                   M为基函数做内积的总刚矩阵，D为求导后的基函数做内积的总刚矩阵
%  输出变量：    U          

%预分配
F=zeros(3*N,1);
%组装右端项
F=M*Q-(D-Q_BR)*R-Q_BL*RL;
%解线性方程
U=M\F;
end