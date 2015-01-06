function [Q_BR,Q_BL] = assemble_Q_Elem(N,h)
%函数 assemble_Q_Elem 组装Q方程里的矩阵
%   输入变量：N为单元数，h为步长
%   输出变量：
%                    Q_BR,Q_BL关于Q方程边界右端点的总刚矩阵,左端点的总刚矩阵

%预分配
bl = zeros(3,3);
br = zeros(3,3);
Q_BR=zeros(3*N,3*N);
Q_BL=zeros(3*N,3*N);

%经过推到得到两个基函数相乘在单元内的右左端点的值只与h有关,无需遍历单元
br=[1,-h/2,h^2/4;h/2,-h^2/4,h^3/8;h^2/4,-h^3/8,h^4/16];
bl=[1,-h/2,h^2/4;-h/2,h^2/4,-h^3/8;h^2/4,-h^3/8,h^4/16];
% 组装
for i=1:N
    Q_BR(3*i-2:3*i,3*i-2:3*i)=br;
    Q_BL(3*i-2:3*i,3*i-2:3*i)=bl;
end

end


