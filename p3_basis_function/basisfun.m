function [phi1,phi2,phi3,phi4] = basisfun(cv1, cv2, v,string)
% v: the points which is mapping from the gauss points 
 mid=(cv1+cv2)/2.0;
 ha=abs(cv2-cv1)/2.0;
switch string
    case 'Mbasis'
        % basis function
       phi1=ones(size(v,1),1);
       phi2=(v-mid)/ha;
       phi3=(3.0*((v-mid)/ha).^2-1.0)/2.0;
       phi4=(5.0*((v-mid)/ha).^3-3.0*((v-mid)/ha))/2.0;
    case 'Dbasis'
        % D basis function
       phi1=zeros(size(v,1),1);
       phi2=ones(size(v,1),1)/ha;
       phi3=(6.0*(v-mid)/(ha*ha))/2.0;
       phi4=(15.0*((v-mid)/ha).^2-3.0)/(2.0*ha);
end




