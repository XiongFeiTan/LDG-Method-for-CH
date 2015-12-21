function [rhs] = initialzieU(Dof,cv1, cv2, p, w)
%% InitialzieU
rhs=zeros(Dof,1);
%% mapping value & jacobi weight
[v, jw]=elemInt(cv1, cv2, p, w);
%% basis function
[phi1,phi2,phi3,phi4]=basisfun(cv1, cv2, v,'Mbasis');
%% initial function
%% example2 accurary test 
u0=0.25.*exp(-abs(v));
%% example3 peakon solution
% u0=zeros(size(v,1),1);
% for i=1:size(v,1)
%      if v(i)>=-20 && v(i)<=10
%        u0(i)=cosh(v(i)+5)/cosh(15);  
%      else
%        u0(i)=cosh(25-v(i))/cosh(15);    
%      end
% end
%% one  element
rhs(1)=sum(phi1.*u0.*jw);
rhs(2)=sum(phi2.*u0.*jw);
rhs(3)=sum(phi3.*u0.*jw);
rhs(4)=sum(phi4.*u0.*jw);
end