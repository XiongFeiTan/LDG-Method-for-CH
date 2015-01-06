function [R,UL] = assemble_R_RH(U,M,R_BR,R_BL,D,N)
%函数 assemble_R_RH 组装R方程里的矩阵
%   输入变量：U初始值，N为单元数，h为步长，D为求导后的基函数做内积的总刚矩阵
%                   R_BR,R_BL关于R方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：R向量，UL为由U变化而来，边值问题处理左端
                  

%预分配
F=zeros(3*N,1);
UL=zeros(3*N,1);

%边界处理
UL(1:3,1)=U(3*N-2:3*N,1);
for i=2:N
    UL(3*i-2:3*i,1)=U(3*i-5:3*i-3,1);
end
%右端项
F=(R_BR-D)*U-R_BL*UL;

%[l,u]=lu(M);
%求解方程
R=M\F;
end
