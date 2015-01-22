%%
%number of element
N=50;
%computing time
time=1;
%number of time steps
NT=N*4;
%time increment
dt=time/NT;
%Interval about the left and right endpoint
a=-25;
b=25;
%meshgrid generation£¬x is node value vector£¬h is element increment£¬Env is matrix(N,2) contains about left and right endpoint value 
[h,Env,x] = MeshGen1D(a,b,N);
%%
%about calculate R,Q equations
[r_br,r_bl,q_br,q_bl,M,D,R_BR,R_BL,Q_BR,Q_BL] = assemble_RQ(Env,N,h);
%calculate U0 initial value
[U] = Initializeu(Env,M,N);
% assemble R equation's right hand item£¬UR is the U to the right
%   translation of a unit vector
[FR,UR] = assemble_R_RH(U,R_BR,R_BL,D,N);
%calculate the R vector (slove a system of linear equations)
R=M\FR;
%%
%about calculate Q equation
%assemble Q equation's right hand item,RL is the R to the left 
%   translation of a unit vector
[FQ,RL] = assemble_Q_RH(U,R,M,Q_BR,Q_BL,D,N);
%calculate the Q vector (slove a system of linear equations)
Q=M\FQ;
%%
%about calculate P equation
[P_BR,P_BL,RU,FP] = assemble_P(U,R,UR,RL,Env,N,h);
%calculate the P vector (slove a system of linear equations)
P=M\FP;
%%
%R and Q composed of the coupling matrix C  assemble
[C]=Coupling_matrix(M,r_bl,r_br,q_bl,q_br,D,N);
%%
%about the QT equation
[QT_BR,QT_BL,RUP,FQT] = assemble_QT(U,R,UR,RL,P,Env,N,h);
[QT]=RK3(M,FQT,dt,Q,C,N,Env,h);
% assign values to Q, for loop
Q=QT;
%%
%time discretization
for tstep=1:NT
    %assemble the coupling system right hand item
    CRH=zeros(6*N,1,'double');
    CRH(1:3*N,1)=M*Q;
    %calculate the U,R vector (slove a system of linear equations)
    X=zeros(6*N,1,'double');
    X=C\CRH;
    U=X(1:3*N,1);
    R=X(3*N+1:6*N,1);
    
   %new UR,RLneed new processing U,R translation of a unit vector
   UR=zeros(3*N,1,'double');
   UR(3*N-2:3*N,1)=U(1:3,1);
    for i=1:N-1
        UR(3*i-2:3*i,1)=U(3*i+1:3*i+3,1);
    end
    
    RL=zeros(3*N,1,'double');
    for i=2:N
    RL(3*i-2:3*i,1)=R(3*i-5:3*i-3,1);
    end
    RL(1:3,1)=R(3*N-2:3*N,1);
   
    %re-sloved P
    [P_BR,P_BL,RU,FP] = assemble_P(U,R,UR,RL,Env,N,h);
    P=M\FP;
    
    %re-sloved QT
    [QT_BR,QT_BL,RUP,FQT] = assemble_QT(U,R,UR,RL,P,Env,N,h);
    [QT]=RK3(M,FQT,dt,Q,C,N,Env,h);
    Q=QT;
end
%%
%calculate U vector 
[U] = Calculate_U(M,Q,D,R,Q_BL,RL,Q_BR,N);
t = linspace (0.0,time,NT+1);
nu=zeros(N,1);
eu=zeros(N,1);

% numerical solution nu
for i=1:N
     if i>=1&&i<=0.4*N
            nu(i,1)=U(3*i-2,1);
    elseif i>0.4*N&&i<0.6*N
             nu(i,1)=0;
    end
end

% exact solution  eu
    for i= 1:N
        if i>=1&&i<=0.4*N
        x(i)=(x(i)+x(i+1))/2;
        eu(i,1) = 0.25*exp(-abs(x(i)-0.25*time));   
        elseif i>0.4*N&&i<0.6*N
       eu(i,1)=0;
        end
    end
% error estimation
error=abs(nu-eu);
L2_error=sqrt(sum((nu-eu).^2)/N)
Lmax_error=max(abs(nu-eu))
% plot 
plot ( x(1:N),nu,'.b:', x(1:N),eu,'+r:');
legend ( 'numrical solution ', 'exact    solution ') ;
grid on ;
u_min=min(min(nu));
u_max=max(max(nu));
axis([a b u_min-0.2,u_max+0.2]);
title ( sprintf ( ' NStep%d,Time%f \n ', N, 1 ) );
xlabel ( '------X------' ) ;
ylabel ( '------U------' );