%%
%单元数
N=160;
%时间
time=1;
%时间步数
NT=400;
% NT=300;
%时间步长
dt=time/NT;
%区间左右端点
a=-25;
b=25;
%产生网格，x为节点值向量，h为单元步长，Env为单元左右端点的值矩阵，大小为(N,2)
[h,Env,x] = MeshGen1D(a,b,N);
%%
%关于求R,Q方程
[r_br,r_bl,q_br,q_bl,M,D,R_BR,R_BL,Q_BR,Q_BL] = assemble_RQ(Env,N,h);
%计算初值U0
[U] = Initializeu(Env,M,N);
%  数值解 nu
% for i=1:N
%     nu(i,1)=U(3*i-2,1);
% end
% 真解 eu
%     for i= 1:N
%         x(i)=(x(i)+x(i+1))/2;
%         eu(i,1) = 0.25*exp(-abs(x(i)));        
%     end
% % 误差估计
% error=abs(nu-eu)
% Lmax_error=max(abs(nu-eu))

% 组装R方程右端项，UR为U往右端平移一个单元的向量
[FR,UR] = assemble_R_RH(U,R_BR,R_BL,D,N);
%求R向量
R=M\FR;
%%
%关于求Q的方程
%组装右端项，RL由R向量向左端平移一个单元的向量
[FQ,RL] = assemble_Q_RH(U,R,M,Q_BR,Q_BL,D,N);
%求Q向量
Q=M\FQ;
%%
%关于求P的方程
[P_BR,P_BL,RU,FP] = assemble_P(U,R,UR,RL,Env,N,h);
%求解线性方程
P=M\FP;
%%
%关于QT
[QT_BR,QT_BL,RUP,FQT] = assemble_QT(U,R,UR,RL,P,Env,N,h);
%欧拉法求解
F_QT=M*Q+dt*FQT;
QT=M\F_QT;
%赋值给Q,循环调用
Q=QT;
%%
%耦合矩阵，组装
[C]=Coupling_matrix(M,r_bl,r_br,q_bl,q_br,D,N);
%%
%时间离散
for tstep=1:NT
    %组装右端项
    CRH=zeros(6*N,1);
    CRH(1:3*N,1)=M*Q;
    %求解U,R
    X=zeros(6*N,1);
    X=C\CRH;
    U=X(1:3*N,1);
    R=X(3*N+1:6*N,1);
    
   %新的UR,RL需新的处理U,R分别向右和左移动一个单元
   UR=zeros(3*N,1);
   UR(3*N-2:3*N,1)=U(1:3,1);
    for i=1:N-1
        UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
    end
    
    RL=zeros(3*N,1);
    for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
    end
    RL(1:3,1)=R(3*N-2:3*N,1);
   
    %重新求P
    [P_BR,P_BL,RU,FP] = assemble_P(U,R,UR,RL,Env,N,h);
    P=M\FP;
    
    %重新求QT
    [QT_BR,QT_BL,RUP,FQT] = assemble_QT(U,R,UR,RL,P,Env,N,h);
    F_QT=M*Q+dt*FQT;
    QT=M\F_QT;
    Q=QT;
end
%%
%利用上面的Q，求解U
[U] = Calculate_U(M,Q,D,R,Q_BL,RL,Q_BR,N);
t = linspace (0.0,time,NT+1);
nu=zeros(N,1);
eu=zeros(N,1);

% 数值解 nu
for i=1:N
    nu(i,1)=U(3*i-2,1);
end

% 真解 eu
    for i= 1:N
        x(i)=(x(i)+x(i+1))/2;
        eu(i,1) = 0.25*exp(-abs(x(i)-0.25*time));        
    end



% 误差估计
error=abs(nu-eu);
L2_error=norm(eu-nu,2);
Lmax_error=max(abs(nu-eu))

% 画图
plot ( x(1:N),nu,'.b:', x(1:N),eu,'+r:');
legend ( 'numrical solution ', 'exact    solution ') ;
grid on ;
u_min=min(min(nu));
u_max=max(max(nu));
axis([a b u_min-0.2,u_max+0.2]);
title ( sprintf ( ' NStep%d,Time%f \n ', N, 1 ) );
xlabel ( '------X------' ) ;
ylabel ( '------U------' );