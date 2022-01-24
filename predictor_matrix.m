function V=predictor_matrix(Xn)

V=zeros(1,size(Xn,1));
V(1,:)=Xn(:,1);   %Re
V(end+1,:)=Xn(:,2); %f+
V(end+1,:)=Xn(:,3); %phi
V=x2fx(V','interaction')'; %mixed terms
V(1,:)=[];
V(end+1,:)=Xn(:,1).^2;   %Re2
V(end+1,:)=Xn(:,2).^2; %f+2
V(end+1,:)=Xn(:,3).^2; %phi2
V(end+1,:)=Xn(:,1).^-0.2;   %Re-0.2
V(end+1,:)=Xn(:,1).^-0.5; %Re-0.5
V(end+1,:)=Xn(:,3).^-1; %phi-1
D=size(V,1);
% names={'Re', 'f', 'phi', 'Re.f', 'Re.phi', 'f.phi','Re2','f2','phi2','Re-0.2','Re-0.5','phi-1'};
