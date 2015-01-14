function [U] = Initializeu(Env,M,N)
%函数 Initializeu 求初值U
%   输入变量：Env单元节点值，N为单元数，M为基函数做内积的总刚矩阵
%   输出变量：U初值

%预分配
U=zeros(3*N,1);
%右端项
F=zeros(3*N,1);

% 单元循环
for i=1:N
    a=Env(i,1);
    b=Env(i,2);
    mid=(a+b)/2;
    
    % 被积函数   
    f1=@(x)(0.25*exp(-abs(x)));
    f2=@(x)(0.25*exp(-abs(x)).*(x-mid));
    f3=@(x)(0.25*exp(-abs(x)).*(x-mid).^2);

    %求积
    quad1=quadrature(f1,a,b);
    quad2=quadrature(f2,a,b);
    quad3=quadrature(f3,a,b);
    %组装右端项
    F(3*i-2:3*i)=[quad1,quad2,quad3]';
end
%[l,u]=lu(M);
%求解线性方程组，求初值
U=M\F;
end


