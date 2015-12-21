function [v, jw]=elemInt(v1, v2, p, w)
% function element Int compute the integral information on a element.
%    v: the points which is mapping from the gauss points 
%    jw: the jacobi*w
a=(v2-v1)/2.0;
b=(v2+v1)/2.0;
v=a*p+b;
jw=w*a;
end