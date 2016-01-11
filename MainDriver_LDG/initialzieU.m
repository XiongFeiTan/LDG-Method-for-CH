function [rhs] = initialzieU(Dof,cv1, cv2, p, w)
%% InitialzieU
rhs=zeros(Dof,1);
%% mapping value & jacobi weight
[v, jw]=elemInt(cv1, cv2, p, w);
% basis function
[phi1,phi2,phi3]=basisfun(cv1, cv2, v,'Mbasis');
%% initial function
%% accurary test example 2
% u0=0.25*exp(-abs(v));
%% example3 peakon solution
u0=zeros(size(v,1),1);
for i=1:size(v,1)
     if v(i)>=-20 && v(i)<=10
       u0(i)=cosh(v(i)+5)/cosh(15);  
     else
       u0(i)=cosh(25-v(i))/cosh(15);    
     end
end
%% two peakon example 4
% u1=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)>=-20 && v(i)<=10
%       u1(i)=2*cosh(v(i)+5)/cosh(15);  
%     else
%       u1(i)=2*cosh(25-v(i))/cosh(15);    
%     end
% end
% u2=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)>=-10 && v(i)<=20
%       u2(i)=cosh(v(i)-5)/cosh(15);  
%     else
%       u2(i)=cosh(35-v(i))/cosh(15);    
%     end
% end
% u0=zeros(size(v,1),1);
% u0=u1+u2;
%% three peakon example 5
% u1=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)>=-20 && v(i)<=10
%       u1(i)=2*cosh(v(i)+5)/cosh(15);  
%     else 
%       u1(i)=2*cosh(25-v(i))/cosh(15);    
%     end
% end
% u2=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)>=-18 && v(i)<=12
%       u2(i)=cosh(v(i)+3)/cosh(15);  
%     else
%       u2(i)=cosh(27-v(i))/cosh(15);    
%     end
% end
% u3=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)>=-16 && v(i)<=14
%       u3(i)=0.8*cosh(v(i)+1)/cosh(15);  
%     else
%       u3(i)=0.8*cosh(29-v(i))/cosh(15);    
%     end
% end
% u0=zeros(size(v,1),1);
% u0=u1+u2+u3;
%% solution with a discontinuous example 6
% u0=10./((3+abs(v)).*(3+abs(v)));
%% breakup of the plateau traveling wave example 7
% u0=zeros(size(v,1),1);
% for i=1:size(v,1)
%     if v(i)<=-5
%       u0(i)=0.6*exp(v(i)+5);  
%     elseif v(i)>=5
%       u0(i)=0.6*exp(-v(i)+5);   
%     else
%       u0(i)=0.6;
%     end
% end
%% one  element
rhs(1)=sum(phi1.*u0.*jw);
rhs(2)=sum(phi2.*u0.*jw);
rhs(3)=sum(phi3.*u0.*jw);
end