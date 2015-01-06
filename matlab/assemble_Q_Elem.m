function [Q_BR,Q_BL] = assemble_Q_Elem(N,h)
%���� assemble_Q_Elem ��װQ������ľ���
%   ���������NΪ��Ԫ����hΪ����
%   ���������
%                    Q_BR,Q_BL����Q���̱߽��Ҷ˵���ܸվ���,��˵���ܸվ���

%Ԥ����
bl = zeros(3,3);
br = zeros(3,3);
Q_BR=zeros(3*N,3*N);
Q_BL=zeros(3*N,3*N);

%�����Ƶ��õ���������������ڵ�Ԫ�ڵ�����˵��ֵֻ��h�й�,���������Ԫ
br=[1,-h/2,h^2/4;h/2,-h^2/4,h^3/8;h^2/4,-h^3/8,h^4/16];
bl=[1,-h/2,h^2/4;-h/2,h^2/4,-h^3/8;h^2/4,-h^3/8,h^4/16];
% ��װ
for i=1:N
    Q_BR(3*i-2:3*i,3*i-2:3*i)=br;
    Q_BL(3*i-2:3*i,3*i-2:3*i)=bl;
end

end

