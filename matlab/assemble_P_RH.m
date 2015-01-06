function [P] = assemble_P_RH(M,P_BR,P_BL,RU,N)
%函数 assemble_P_Elem  求P向量
%   输入变量：N为单元数,M为基函数做内积的总刚矩阵
%                    P_BR,P_BL关于P方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：P向量
%                   

%预分配
F=zeros(3*N,1);
%右端项
F=-RU+P_BR-P_BL;
%[l,u]=lu(M);
%求解线性方程
P=M\F;
end