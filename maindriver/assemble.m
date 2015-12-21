function [C,M]= assemble(Dof,N,m,d,jr,j,jl )
M=kron(speye(N),m);
C1=kron(sparse(diag(ones(N-1,1),-1)),jl)+kron(sparse(diag(ones(N-1,1),1)),jr)+kron(speye(N),(d+j));
C1(1:Dof,Dof*(N-1)+1:Dof*N)=jl;
C1(Dof*(N-1)+1:Dof*N,1:Dof)=jr;
C=[M,C1;C1,M];
end

