function [U,R]=Coupling_matrix(Q,M,Q_BL,Q_BR,R_BL,R_BR,D,N)
%coupling matrix equation
%calculate U,R
%predistribution
C=zeros(6*N,6*N);
c1=zeros(3*N,3*N);
c2=zeros(3*N,3*N);
c3=zeros(3*N,3*N);
F1=zeros(3*N,1);
F2=zeros(3*N,1);
X=zeros(6*N,1);

%
c1=M;
%c2
for i=1:N-1
    c2(3*i-2:3*i,3*i+1:3*i+3)=-Q_BR(3*i-2:3*i,3*i-2:3*i);  
    c2(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)+Q_BL(3*i-2:3*i,3*i-2:3*i);
    c3(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)-Q_BR(3*i-2:3*i,3*i-2:3*i);
end
c2(3*N-2:3*N,1:3)=-Q_BR(3*N-2:3*N,3*N-2:3*N);
c2(3*i-2:3*i,3*i-2:3*i)=D(3*i-2:3*i,3*i-2:3*i)+Q_BL(3*i-2:3*i,3*i-2:3*i);
c2(3*N-2:3*N,3*N-2:3*N)=D(3*N-2:3*N,3*N-2:3*N)+Q_BL(3*N-2:3*N,3*N-2:3*N);
c3(3*N-2:3*N,3*N-2:3*N)=D(3*N-2:3*N,3*N-2:3*N)-Q_BR(3*N-2:3*N,3*N-2:3*N);
%
for i=2:N
    c3(3*i-2:3*i,3*i-5:3*i-3)=R_BL(3*i-2:3*i,3*i-2:3*i);
end
c3(1:3,3*N-2:3*N)=R_BL(3*N-2:3*N,3*N-2:3*N);

%
C(1:3*N,1:3*N)=c1;
C(3*N+1:6*N,3*N+1:6*N)=c1;
C(1:3*N,3*N+1:6*N)=c2;
C(3*N+1:6*N,1:3*N)=c3;
%

F1=M*Q;
F(1:3*N,1)=F1;
F(3*N+1:6*N,1)=F2;

%
X=C\F;
U=X(1:3*N,1);
R=X(3*N+1:6*N,1);

end

