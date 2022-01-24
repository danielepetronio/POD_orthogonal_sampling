function [U_red V_red A_mean]=POD(A,threshold)

%% POD
A_mean=mean(A,2); %%mean cpt distribution
A_f=A-A_mean; %% cpt fluctuations
C=A_f'*A_f; %% cross-correlation
[V,eps2]=eig(C);
eps2=diag(eps2)'; eps2=fliplr(eps2); %% eigenvalues (or squared svd singular values)
V=fliplr(V); %% eigenvectors (POD coefficients)
U=A_f*V; %% moded
Etot=sum(eps2); %%total energy of the dataset
E=cumsum(eps2)/Etot; %%cumulative energy

%% model-reduction (99% energy)
ind_in=find(E<threshold); 
d_red=length(ind_in); %% reduced order (approximated rank)
U_red=U(:,ind_in); %%reduced modes
V_red=V(:,ind_in); %%reduced coefficients