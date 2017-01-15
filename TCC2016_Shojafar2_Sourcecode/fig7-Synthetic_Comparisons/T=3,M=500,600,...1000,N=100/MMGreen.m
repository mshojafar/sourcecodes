% total_energyt1=xlsread('Lyap_total_energy_PMR1.5.csv');
%     csvwrite('Lyap_total_energy_PMR1.5.dat',total_energyt1);
%     
%     tot_reconfig_Energy1=xlsread('Lyap_tot_reconfig_Energy_PMR1.5.csv');
%     csvwrite('Lyap_tot_reconfig_Energy_PMR1.5.dat',tot_reconfig_Energy1);
%     
%     Opt_Comm_Energyt1=xlsread('Lyap_tot_Comm_Energy_PMR1.5.csv');
%     csvwrite('Lyap_tot_Comm_Energy_PMR1.5.dat',Opt_Comm_Energyt1);
%     
%     Opt_Comp_Energyt1=xlsread('Lyap_tot_Comp_Energy_PMR1.5.csv');
%     csvwrite('Lyap_tot_Comp_Energy_PMR1.5.dat',Opt_Comp_Energyt1);
%     
%     save('Lyap_ScnerioM50_5_100N100F1T3.mat','total_energyt1','tot_reconfig_Energy1',...
%         'Opt_Comm_Energyt1','Opt_Comp_Energyt1');

clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%difinite variables
M=1000;% numbrs of virtual machines
M_vector=[500:100:1000];
N=100;% workload numbers
Q=3;% numbers of Discrete ranges for Frequency of each VM
%T=0.1*ones(1,M);% seconds/ time for each Vm for computation
T=3;
T_t=5;
%T_vector=[0.1:0.1:0.9].*T_t;
%T_vector_range=[0.1:0.1:0.9];
%T_vector=[1:0.1:T_t-1];
R_tot=10000;% Mbit / seconds
temp=0;
ALPHA=[0.1,0.5,0.9];

k_e=[0.005]*1000;% Joule/(Mega Hz)^2 % scenario 1
%Discrete_F=[0.15, 1.867, 2.133,2.533,2.668];%(GhZ)% for Interl Nehalem Quad-core Processor in INFOCOM 2009 G.V. Laszewski et al
% Processing time	0.08 (s)
% Raw Image Size	24000000.00  (bits)
% Compression	23.00
% Compressed Image size	1043478.26 (bits)
% CPU [0.3, 0.533, 0.667, 0.800, 0.933] %(GhZ)%==[1.3, 2.32, 2.94, 3.48, 4.06] %(Mbit/s)% Processing Rate
% CPU [0.15, 1.867, 2.133,2.533,2.668] %(GhZ)%==[0.65217, 8.12, 9.27, 11.01, 11.60] %(Mbit/s)% Processing Rate

%Discrete_F=[0.65217, 8.12, 9.27, 11.01, 11.60];%(Mbit/s)% for Interl Nehalem Quad-core Processor in INFOCOM 2009 G.V. Laszewski et al
%Discrete_F=[1.3, 2.32, 2.94, 3.48, 4.06];%(Mbit/s)% for power-scalable real cluster Crusoe with CPU type TM-5800 in CLUSTER 2006 in http://dx.doi.org/10.1109/CLUSTR.2006.311839
Discrete_F=[8.12, 9.27, 11.01, 11.60];%(Mbit/s)% for Interl Nehalem Quad-core Processor in INFOCOM 2009 G.V. Laszewski et al
omega=1;% cost portion of VM consumption
%P_Idle=0.5*ones(M,1);% Watt: Idle Power for each VM
%E_max=100*ones(1,M);% Joule-maximum tolarated energy for each VM
%Delta=0.1*ones(1,M); %seconds-considered time of processing for each VM
W=1;% MHZ frequency needed for communication channel
Omega=0.005;%milli Watt


%L_tot=60 + (80-60).*rand(1,1000); %carico uniformemente distribuito   70+-10 with avg(L_{tot})=70 (PMR=1.14285) %UTILIZZATO PER Scenario 2
%save('workload1000M.mat','L_tot');
% L_tot_Vector = load('workload1000M','L_tot');
% L_tot=L_tot_Vector.L_tot;
%L_tot_Vector  = load('workload1000');
%L_tot=L_tot_Vector.workload1000;
% L_tot_Vector  = load('UNI1_IMC2010_DC1');
% L_tot_Vector  = load('UNI1_IMC2010_DC2');
% L_tot_Vector  = load('UNI1_IMC2010_DC3');
% L_tot_Vector  = load('UNI1_IMC2010_DC4');
% L_tot_Vector  = load('UNI2_IMC2010_DC1');
% L_tot_Vector  = load('UNI2_IMC2010_DC2');
% L_tot_Vector  = load('UNI2_IMC2010_DC3');
% L_tot_Vector  = load('UNI2_IMC2010_DC4');
%L_tot=L_tot_Vector.L_tot(1:N);
%pmr=(max(L_tot)-mean(L_tot))/mean(L_tot);
%% synthetic workload definitions
mL_tot=9000;
pmr=1.5;
a=mL_tot*(pmr-1);
if pmr>1
%L_tot=(mL_tot-a) + 2.*(a).*rand(1,N);
%save('Ltot900PMR1.5.mat','L_tot');
%save('Ltot9000PMR1.5.mat','L_tot');

%L_tot=load('Ltot900PMR1.5.mat','L_tot');

%save('Ltot90000PMR1.5.mat','L_tot');
%L_tot=load('Ltot90000PMR1.5.mat','L_tot');
L_tot=load('Ltot9000PMR1.5.mat','L_tot');


L_tot=L_tot.L_tot(1:N);
else
    L_tot=(mL_tot-a).*ones(1,N);
end
for i=1:N
    if (L_tot(i)>0)
        L_tot(i)=L_tot(i);
    else
        L_tot(i)=0;
    end
end
%L_tot=1400*ones(1,N);
%M=ceil(max(L_tot)./(T(1)*max(Discrete_F)))+1;
%L_tot=L_tot./T_t;
%L_tot=M*max(Discrete_F)/max(L_tot).*L_tot(1:N); % to make L_tot tune to work with
%high range frequency
%L_tot = 8 + PMR .* rand(1,N);%task-size of each VM

f_opt_wokloads=zeros(N,M); % optimal frequencies of M VMs for each workload
%L_opt_wokloads=zeros(N,M); % optimal Load of M VMs for each workload
t_opt_wokloads=zeros(N,Q+1); % optimal rate of M VMs for each workload

%%%%%%% CVX solution
%Opt_reconfig_Energy1=zeros(length(L_tot),1);
MM=M_vector;
Opt_Comp_Energy=zeros(length(L_tot),length(MM));
Opt_Comp_Energy1=zeros(length(L_tot),length(MM));

Opt_reconfig_Energy=zeros(length(L_tot),length(MM));
Opt_reconfig_Energy11=zeros(length(L_tot),length(MM));
Opt_reconfig_Energy12=zeros(length(L_tot),length(MM));
Opt_reconfig_Energy13=zeros(length(L_tot),length(MM));

Opt_Comm_Energy=zeros(length(L_tot),length(MM));
Opt_Comm_Energy1=zeros(length(L_tot),length(MM));
Opt_Comm_Energy2=zeros(length(L_tot),length(MM));
Opt_Comm_Energy3=zeros(length(L_tot),length(MM));
Opt_Comm_Energy4=zeros(length(L_tot),length(MM));

total_energy=zeros(length(L_tot),length(MM));
total_energy1=zeros(length(L_tot),length(MM));
total_energy2=zeros(length(L_tot),length(MM));
total_energy3=zeros(length(L_tot),length(MM));
total_energy4=zeros(length(L_tot),length(MM));


tic

 %Total_Time=zeros(length(MM),length(L_tot),M,Q+1);
 t_opt_wokloads_per_Q=zeros(1,Q+1);
 SLA_vioaltion=zeros(length(MM),length(L_tot));
 %%
    %%code body
     Opt_reconfig_Energy1=zeros(length(L_tot),length(MM));
     
     
  for jj=1:length(MM),
    %T=T_vector(jj);
    
    j=M_vector(jj);
    temp=0;
    f_zero=0;
    F=ones(j,1)*Discrete_F; % Discrete-frequency vector for M vms (Hemogenous VMs)
    PowerFunc=pow_p(F,3);
    P_Idle=0.005*ones(j,1);% Watt: Idle Power for each VM channel
    f_max=max(Discrete_F).*ones(1,j);% maximum frequencies for each VM
    Delta_F=zeros(j,1);
    
    Total_Time=zeros(length(MM),length(L_tot),j,Q+1);
    
    
    for i=1:(length(L_tot)),
        tic
       if L_tot(i)> sum(f_max(:).*T)
           SLA_vioaltion(jj,i)=1;
           
           Opt_Comp_Energy(i,jj)=NaN;
           Opt_Comm_Energy1(i,jj)=NaN;
           Opt_reconfig_Energy11(i,jj)=NaN;
           total_energy1(i,jj)=NaN;
           
           Opt_Comp_Energy(any(isnan(Opt_Comp_Energy),2),:) = [];
           Opt_Comm_Energy1(any(isnan(Opt_Comm_Energy1),2),:) = [];
           Opt_reconfig_Energy11(any(isnan(Opt_reconfig_Energy11),2),:) = [];
           total_energy1(any(isnan(total_energy1),2),:) = [];
          

           continue
          %error('Problem SLA violation')
       end
       
       if L_tot(i)==0
           Total_Time( ~any(Total_Time,2), : ) = []; % delete zero rows in matlab.
          continue
       end
                       
        %*******************Feasibility condition************************************

        condizione_feas_workloadcomputation= (L_tot(i)<= sum(f_max(:).*T));
        condizione_feas_channel=(L_tot(i)<= R_tot*Q*(T_t-T/2));
        
        feasibility=condizione_feas_workloadcomputation && condizione_feas_channel;
        
        if ~feasibility
            'VM='
            M
            error('Problem unfeasible')
        else
            %  'Problem Feasibile!'
        end
    %**************************************************************************
    cvx_begin quiet
    cvx_solver sedumi
    variable t_opt(j,Q+1) nonnegative; % optimum load for each VM
    RTT=70; %Ms
    minimize ((omega(1)*sum(PowerFunc(:).*t_opt(:)))+(sum(Omega*(RTT.*2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T))).^2)+sum(P_Idle(:))));
    
    %s.t. constrainsts
    subject to
    sum(sum(F(:,2:Q+1).*t_opt(:,2:Q+1)))-L_tot(i)==0;
    for  kk=1:j
    0<=sum(t_opt(kk,:))-T<=0;
   % (2*(T_t-T(1))).*(sum(F(kk,:).*t_opt(kk,:))./(sum(F(kk,:).*t_opt(kk,:),2)))+T(1)-T_t<=0;
    end
   
    0<= t_opt <= T;

    0<= F <= max(Discrete_F);

    sum((sum(F(:,1:Q+1).*t_opt(:,1:Q+1),2))/(T_t-T))-R_tot<=0;
    
   cvx_end
    %Reconfiguration cost******begin

    if (t_opt(:,1)>0)
    Delta_F(:)=Delta_F(:)+((F(:,1)-f_zero(:)).^2); % for saving the first Delta-time for each VM 
    else
    Delta_F=zeros(j,1); % for saving the first Delta-time for each VM if the first t_opt for each VM is zero: VM no in idle mode at starting finding t_opt
    end
    for  l=2:(Q+1)
        if (t_opt(:,l)>0)
        Delta_F(:)=Delta_F(:)+((F(:,l)-F(:,l-1)).^2);
        f_zero=F(:,l); % for saving the last active Frequency for each VM for the next incoming workload
        end
    end
    Opt_reconfig_Energy1(i,jj)=sum(Delta_F(:));
    Delta_F=zeros(j,1);
    %%55
    %%Results
    t_opt_wokloads(i,:)=t_opt(j,:);
    t_opt_wokloads(any(isnan(t_opt_wokloads),2),:) = [];
    Total_Time(jj,i,:,:)=t_opt;
    if isnan(t_opt) 
      %  Total_Time(any(isnan(Total_Time),2),:) = [];
     t_opt= Total_Time(jj-1,i,:,:);
        continue
    end
    
    if isnan(t_opt) 
        Total_Time(any(isnan(Total_Time),2),:) = [];
        continue
    end
    
    
     Opt_Comp_Energy(i,jj)=omega(1)*sum(PowerFunc(:).*t_opt(:));
      
     Opt_reconfig_Energy11(i,jj)=k_e(1)*Opt_reconfig_Energy1(i,jj);
  
     Opt_Comm_Energy1(i,jj)=(sum(Omega*(RTT.*2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T))).^2)+sum(P_Idle(:)));    
     %Opt_Comm_Energy1(i,jj)=(T_t-T(1))*((Zeta5*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
    
     total_energy1(i,jj)=Opt_Comp_Energy(i,jj)+Opt_reconfig_Energy11(i,jj)+Opt_Comm_Energy1(i,jj);
  
    
     %cvx_cputime_1(i,j)=toc;
    cvx_cputime_2(i,jj)=cvx_cputime;
    %f_opt_wokloads(i,:)=f(1);
    %L_opt_wokloads(i,:)=L(1);
    
    % t_opt_wokloads(~any(isnan(t_opt_wokloads)));
    % t_opt_wokloads(isnan(t_opt_wokloads))=0;
    
    end % for jobs (workloads)
  end % for VMs
total_energyt1=zeros(1,length(MM));

Opt_Comp_Energyt1=zeros(1,length(MM));

Opt_Comm_Energyt1=zeros(1,length(MM));

tot_reconfig_Energy1=zeros(1,length(MM));

time_meanWL=zeros(length(MM),M,Q+1);
time_meanWL_M=zeros(length(MM),Q+1);
sum_time_meanWL_M=zeros(1,length(MM));

% for jj=1:length(MM)
%     for i=1:length(L_tot)
%         if (L_tot(i)==0)
%             Total_Time( ~any(Total_Time,2), : ) = []; % delete zero rows in matlab.
%             for j=1:M
%                 for qq=1:Q+1
%                     
%                 end
%             end
%         end
%         
%     end
% end
%Total_Time( ~any(Total_Time,2), : ) = []; % delete zero rows in matlab.

for jj=1:length(MM)
    Opt_Comp_Energyt1(jj)=mean(Opt_Comp_Energy(:,jj));
    Opt_Comm_Energyt1(jj)=mean(Opt_Comm_Energy1(:,jj));
    tot_reconfig_Energy1(jj)=mean(Opt_reconfig_Energy11(:,jj));
    total_energyt1(jj)=mean(total_energy1(:,jj));
end
   

for jj=1:length(MM)
    for j=1:M
        for qq=1:Q+1
            time_meanWL(jj,j,qq)=mean(Total_Time(jj,:,j,qq));
        end
    end
end
for jj=1:length(MM)
    for  qq=1:Q+1
        time_meanWL_M(jj,qq)=mean(time_meanWL(jj,:,qq));
    end
end
for jj=1:length(MM)
    sum_time_meanWL_M(jj)=sum(time_meanWL_M(jj,:));
end
save('MMGreen_ScnerioM500_100_1000N100F1T2.mat','total_energyt1','tot_reconfig_Energy1',...
        'Opt_Comm_Energyt1','Opt_Comp_Energyt1','Total_Time','time_meanWL_M', 'time_meanWL');

    xlswrite('total_energy_PMR1.5.xlsx',total_energyt1);
    csvwrite('total_energy_PMR1.5.dat',total_energyt1);
    
    xlswrite('tot_reconfig_Energy_PMR1.5.xlsx',tot_reconfig_Energy1);
    csvwrite('tot_reconfig_Energy_PMR1.5.dat',tot_reconfig_Energy1);
    
    xlswrite('tot_Comm_Energy_PMR1.5.xlsx',Opt_Comm_Energyt1);
    csvwrite('tot_Comm_Energy_PMR1.5.dat',Opt_Comm_Energyt1);
    
    xlswrite('tot_Comp_Energy_PMR1.5.xlsx',Opt_Comp_Energyt1);
    csvwrite('tot_Comp_Energy_PMR1.5.dat',Opt_Comp_Energyt1);
    
    
xlswrite('time_T_T_tperQ.xlsx',time_meanWL_M);
csvwrite('time_T_T_tperQ.dat',time_meanWL_M);

%xlswrite('Time.xlsx',Total_Time);
%csvwrite('csvtime.dat',Total_Time);
time=toc;

Results=load('MMGreen_ScnerioM500_100_1000N100F1T2.mat','total_energyt1','tot_reconfig_Energy1',...
        'Opt_Comm_Energyt1','Opt_Comp_Energyt1','Total_Time','time_meanWL_M', 'time_meanWL');    
total_energyt1=Results.total_energyt1;
Opt_Comm_Energyt1=Results.Opt_Comm_Energyt1;
Opt_Comp_Energyt1=Results.Opt_Comp_Energyt1;
tot_reconfig_Energy1=Results.tot_reconfig_Energy1;
T1=Results.time_meanWL;
T2=Results.time_meanWL_M;


%%
%%%%%%% simulation plots

figure(2000) %fig 7a
plot(M_vector,total_energyt1,'-bs','markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}_{tot}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\overline{L}_{tot},\:PMR=1.5$','Interpreter','latex','FontSize',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);
grid on

figure(2001) %fig 7b
plot(M_vector,tot_reconfig_Energy1,'-bs','markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}_{REC}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\overline{L}_{tot},\:PMR=1.5$','Interpreter','latex','FontSize',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);

grid on

figure(2002) %fig 7c
plot(M_vector,Opt_Comp_Energyt1,'-bs','markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}_{CPc}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\overline{L}_{tot},\:PMR=1.5$','Interpreter','latex','FontSize',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);
grid on

figure(2003) %fig 7d
plot(M_vector,Opt_Comm_Energyt1,'-bs','markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}^{CMc}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\overline{L}_{tot},\:PMR=1.5$','Interpreter','latex','FontSize',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);
grid on

figure(2003) %fig 5a
plot(M_vector,time_meanWL_M(:,1),'-bs',M_vector,time_meanWL_M(:,2),'-kd',...
     M_vector,time_meanWL_M(:,3),'-ro',M_vector,time_meanWL_M(:,4),'-k*',...
     'markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$Agreegate\:Time\:(\%)$','Interpreter','latex','FontSize',40);
legend('$F_{i1}$','$F_{i2}$','$F_{i3}$','$F_{iQ}$',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);
grid on


SLA_vioaltiont1=zeros(1,length(T_vector_range));

for i=1:length(T_vector_range)
SLA_vioaltiont1(i)=sum(SLA_vioaltion(1,i,:));
end
figure(2004) %fig 5g- SLA violation over PMRs
plot(M_vector,SLA_vioaltiont1,'-bs','markers',20,'LineWidth',6);
xlabel('$M$','Interpreter','latex','FontSize',40);
ylabel('$Aggregate SLA Vilotation$','Interpreter','latex','FontSize',40);
legend('$PMR=1.5$','Interpreter','latex','FontSize',40);
title(['M= ' num2str(M) ', N= ' num2str(N) ', F_Q= ' num2str(max(Discrete_F))],'Interpreter','latex','FontSize',30);
grid on

