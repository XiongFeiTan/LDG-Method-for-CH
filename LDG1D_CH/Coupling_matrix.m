function [C]=Coupling_matrix(M,r_bl,r_br,q_bl,q_br,D,N)
%函数 Coupling_matrix 组装耦合方程里的矩阵
%   输入变量：N为单元数
%   输出变量：C为耦合矩阵

%预分配                  
C=zeros(6*N,6*N);
%c1=zeros(3*N,3*N);
c2=zeros(3*N,3*N);
c3=zeros(3*N,3*N);

%c3
for i=1:N-1
    c3(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)+r_bl;
    c3(3*i-2:3*i,3*i+1:3*i+3)=-r_br;
end
c3(3*N-2:3*N,3*N-2:3*N)=D(3*N-2:3*N,3*N-2:3*N)+r_bl;
c3(3*N-2:3*N,1:3)=-r_br;

%c2
for i=2:N
    c2(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)-q_br;
    c2(3*i-2:3*i,3*i-5:3*i-3)=q_bl;
end
c2(1:3,1:3)=D(1:3,1:3)-q_br;
c2(1:3,3*N-2:3*N)=q_bl;

%组装耦合矩阵
C(1:3*N,1:3*N)=M;
C(3*N+1:6*N,3*N+1:6*N)=M;
C(1:3*N,3*N+1:6*N)=c2;
C(3*N+1:6*N,1:3*N)=c3;
end


