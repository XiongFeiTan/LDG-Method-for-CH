function [Q3] = RK3(tDof,Q,QT,dt,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4)
  % 3rd order SSP Runge-Kutta
  %%
  % SSP RK Stage 1.
  Q1  = Q  + dt.*QT;
  %%
  % SSP RK Stage 2.
  [QT2]  = getQT(tDof,Q1,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4);
  Q2   = (3*Q  + Q1  + dt.*QT2 )/4;
%%
  % SSP RK Stage 3.
  [QT3]  = getQT(tDof,Q2,M,C,m,N,jw,phi1,phi2,phi3,phi4,dphi1,dphi2,dphi3,dphi4);
  Q3  = (Q  + 2*Q2  + 2*dt.*QT3 )/3;
end
