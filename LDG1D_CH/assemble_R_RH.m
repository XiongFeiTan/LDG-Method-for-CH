function [FR,UR] = assemble_R_RH(U,R_BR,R_BL,D,N)
%函数 assemble_R_RH 组装R方程里的矩阵
%   输入变量：U初始值，N为单元数，h为步长，D为求导后的基函数做内积的总刚矩阵
%                   R_BR,R_BL关于R方程边界右端点的总刚矩阵,左端点的总刚矩阵
%   输出变量：FR向量，UL为由U变化而来，边值问题处理左端
                 
%预分配
FR=zeros(3*N,1);
UR=zeros(3*N,1);
%边界处理
UR(3*N-2:3*N,1)=U(1:3,1);
for i=1:N-1
    UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
end
%右端项
FR=-(R_BL+D)*U+R_BR*UR;
end
