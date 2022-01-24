clc;
clear all;
close all;

load('all_data.mat') %%load data
N=size(X,2); %%number total tests
omega_meas=trapz(y,CPT,1)*100; %%measured omega(profile-losses)

%% starting point
ind_in=[3,9,28,43]; %%corner of the text matrix
start=length(ind_in); %%size of the initial point

rm=0; %%select 0 for regression method: Gaussian Process
% rm=1; %%select 1 for regression method: Lasso

fig1=figure; %% figure test matrix filling
x0=50; y0=50; width=1700; height=600;
set(gcf,'position',[x0,y0,width,height])
ax=subplot(1,3,1);
plot3(X(3,:),X(1,:),X(2,:),'bo','MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor',[0,0,1]); hold on;
set(ax,'XGrid','on'); set(ax,'YGrid','on'); set(ax,'ZGrid','on'); set(ax,'View',[-20 30]);
set(ax,'FontName','Times New Roman','FontSize',19);
set(ax,'plotboxaspectratio',[1,1.113,1.03],'Box','on');
xlabel('φ','FontName','Times New Roman','FontSize',26,'FontAngle','Italic','verticalalignment','middle')
ylabel('Re','FontName','Times New Roman','FontSize',26,'FontAngle','Italic','verticalalignment','bottom')
zlabel('f^+','FontName','Times New Roman','FontSize',26,'FontAngle','Italic','rotation',0)


%% (predictive analysis)
for k=start:N    
    x_in=X(:,ind_in); %%in data at iteration k 
    
    %%%upload figure test matrix
    ax=subplot(1,3,1); hold on;
    plot3(X(3,ind_in),X(1,ind_in),X(2,ind_in),'o','MarkerSize',11,'MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0]);

    %%POD
    [U_red V_red cpt_mean]=POD(CPT(:,ind_in),0.99); %%POD with the only in tests
    xref=sum(x_in.^2,2); x_in_norm=x_in./sqrt(xref); Xn=X./sqrt(xref); %% normalization inputs for the fitting
    e = {'eig1','eig2','eig3','eig4','eig5','eig6','eig7','eig8','eig9','eig10'}; e = e(1:size(U_red,2));

%%%%%%%%%%% regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rm==0 %%% Gaussian Process regression %%%
    for i=1:size(U_red,2)
        yI=-V_red(:,i);
        gprMdl =fitrgp(x_in_norm',yI,'KernelFunction','ardsquaredexponential');
        GPR.(e{i}) = gprMdl; clear yI gprMdl;
    end
    % prediction over the entire test matrix (gaussian process)
    for i=1:size(U_red,2)
        [Vmod1 Vmod_std1] = predict(GPR.(e{i}),Xn');
        Vmod(:,i)=Vmod1; Vmod_std(:,i)=Vmod_std1; clear Vmod1 Vmod_std1;
    end
    elseif rm==1 %%% Lasso regression %%%
    Vn=predictor_matrix(x_in_norm');
    for i=1:size(U_red,2)
        yI=-V_red(:,i);
        [alfa,FitInfo] = lasso(Vn',yI,'Alpha',1,'Intercept',true,'CV',k);
        idx=FitInfo.Index1SE; %%choose min SE of k-fold (leave-one-out) cross-validation
        alfa_0(i)=FitInfo.Intercept(idx);
        alfa_tot(:,i)=alfa(:,idx);
        fi_in=Vn'*alfa_tot(:,i)+alfa_0(i);
        sigma(i)=sqrt(sum((yI-fi_in).^2)/k);
        V1=Vn; V1(find(alfa_tot(:,i)==0),:)=[]; V1(end+1,:)=1;
        covariance.(e{i}) = sigma(i)^2*inv(V1*V1');
        clear alfa YI FitInfo V1 fi_in
    end
    % prediction over the entire test matrix (lasso)
    Vn=predictor_matrix(Xn');
    Vmod=Vn'*alfa_tot+alfa_0;  %%add covariance for estimating Vmod_std
    for i=1:size(U_red,2)
        V1=Vn; V1(find(alfa_tot(:,i)==0),:)=[]; V1(end+1,:)=1;
        Vmod_std(:,i)=diag(sqrt(V1'*covariance.(e{i})*V1+sigma(i)^2)); clear V1;
    end
    clear sigma covariance;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cpt_fit=-U_red*Vmod'+cpt_mean;  %% fitted cpt distributions over the entire test matrix
    omega_fit=trapz(y,cpt_fit,1)*100; %% fitted omega (profile-losses)
    for i=1:size(U_red,2) %%calculation of std of omega
        sigma(:,:,i)=-U_red(:,i)*Vmod_std(:,i)';
    end
    cpt_std=sqrt(sum(sigma.^2,3)); clear sigma; %%std over cpt distribution
    deltay=y(2)-y(1);
    omega_std=sqrt(sum(cpt_std.^2,1))*deltay*100; %%std over omega
    
    %%errors
    error=abs(omega_fit-omega_meas);
    MEANerror(k)=mean(error./omega_meas)*100;
    MAXerror(k)=max(error./omega_meas)*100;

    %%standard deviations
    MEANstd(k)=mean(omega_std./omega_meas)*100;
    MAXstd(k)=max(omega_std./omega_meas)*100;
    
    %%error and std plots(figures 10 nad11 in the paper)
    ax=subplot(2,3,2); hold on;
    if rm==0
        Color=[0.1 0.8 0.1];
    elseif rm==1
        Color=[0.8 0.1 0.1];
    end
    plot(MAXerror,'-s','MarkerEdgeColor',Color,'MarkerFaceColor',Color,'Color',Color,'LineWidth',1.5,'Markersize',4.5); xlim([4 N]);
    ylim([0 30]); set(ax,'XGrid','on','GridColor','k'); set(ax,'YGrid','on','GridColor','k');
    set(ax,'FontName','Times New Roman','fontsize',16);
    set(ax,'Color',[1 1 1],'Box','on');
    xlabel('n','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')
    ylabel('E_m_a_x(ω)','FontName','Times New Roman','FontSize',22,'FontAngle','Italic','rotation',90)
    ax=subplot(2,3,3); hold on;
    plot(MEANerror,'-s','MarkerEdgeColor',Color,'MarkerFaceColor',Color,'Color',Color,'LineWidth',1.5,'Markersize',4.5); xlim([4 N]);
    ylim([0 12]); set(ax,'XGrid','on','GridColor','k'); set(ax,'YGrid','on','GridColor','k');
    set(ax,'FontName','Times New Roman','fontsize',16);
    set(ax,'Color',[1 1 1],'Box','on');
    xlabel('n','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')
    ylabel('E_m_e_a_n(ω)','FontName','Times New Roman','FontSize',22,'FontAngle','Italic','rotation',90)
    ax=subplot(2,3,5); hold on;
    plot(MAXstd,'-s','MarkerEdgeColor',Color,'MarkerFaceColor',Color,'Color',Color,'LineWidth',1.5,'Markersize',4.5); xlim([4 N]);
    ylim([0 5]); set(ax,'XGrid','on','GridColor','k'); set(ax,'YGrid','on','GridColor','k');
    set(ax,'FontName','Times New Roman','fontsize',16);
    set(ax,'Color',[1 1 1],'Box','on');
    xlabel('n','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')
    ylabel('σ_m_a_x(ω)','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')
    ax=subplot(2,3,6); hold on;
    plot(MEANstd,'-s','MarkerEdgeColor',Color,'MarkerFaceColor',Color,'Color',Color,'LineWidth',1.5,'Markersize',4.5); xlim([4 N]);
    ylim([0 5]); set(ax,'XGrid','on','GridColor','k'); set(ax,'YGrid','on','GridColor','k');
    set(ax,'FontName','Times New Roman','fontsize',16);
    set(ax,'Color',[1 1 1],'Box','on');
    xlabel('n','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')
    ylabel('σ_m_e_a_n(ω)','FontName','Times New Roman','FontSize',22,'FontAngle','Italic')

    %% identification of the next test (most orthogonal in the POD space, see also algorithm 1 in the paper)
    Vmod_norm=Vmod./sqrt(sum(Vmod.^2,2)); 
    Z=abs(Vmod_norm*Vmod_norm'); %%cross-correlation with modelled projected data
    s=sum(Z(ind_in,:),1); s(ind_in)=100; %%find min correlation (most orthogonal)
    ind_next=find(s==min(s)); ind_next=ind_next(1); %% index of next test   
    ind_in=[ind_in ind_next]; %% upload the cicle

    if k==15 %%figure modelled cpt with a reduced number of tests (figure 12 in the paper)
        fig12=figure_cpt(y,cpt_fit,CPT);
        error_example=error([34,33,28,7])./omega_meas([34,33,28,7])*100;
        std_example=omega_std([34,33,28,7])./omega_meas([34,33,28,7])*100;
        saveas(fig12,'example_cpt.jpg'); close(fig12)
    end

%%%%% create gif for the test matrix filling %%%%%%%%%%%%%%%%%%%%%%
    drawnow 
    frame = getframe(fig1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k==start
        imwrite(imind,cm,'test_matrix_filling.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'test_matrix_filling.gif','gif','WriteMode','append');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear Z s ind_next Vmod_norm    
    clear error omega_fit omega_std cpt_fit cpt_std GPR e x_in x_in_norm x_out xref Vmod Vmod_std alfa_tot alfa_0
end



