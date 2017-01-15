clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%difinite variables
M=10;% numbrs of virtual machines
N=1000% workload numbers
Q=4% numbers of Discrete ranges for Frequency of each VM
T=5*ones(1,M);% seconds/ time for each Vm for computation
%T_vector=[4,5,6];
R_tot=100;% Mega bite / seconds
%R_tot_vector=[10,100];
PMR=1.25;
T_t=7;
temp=0;
ALPHA=[0.1,0.5,0.9];

k_e=[0.005;0.05];% Joule/(Mega Hz)^2 % scenario 1
%k_e=[0.005*1000;0.05*1000;0.5*1000];% Joule/(Giga Hz)^2 % scenario 2
Discrete_F2=[0.15, 1.867, 2.133,2.533,2.668];% for Interl Nehalem Quad-core Processor in INFOCOM 2009 G.V. Laszewski et al
Discrete_F1=[0.3, 0.533, 0.667, 0.800, 0.933];% for power-scalable real cluster Crusoe with CPU type TM-5800 in CLUSTER 2006 in http://dx.doi.org/10.1109/CLUSTR.2006.311839
%F=ones(M,1)*Discrete_F; % Discrete-frequency vector for M vms (Hemogenous VMs)

%F1=ones(M,1)*Discrete_F(2:length(Discrete_F)); % Discrete-frequency vector for M vms (Hemogenous VMs) without F_0 which is belongs to Idle state
%Delta_F=zeros(M,1);  % (f(new)-f(old))^2
%f_zero=zeros(M,1);% starting frequncy of each VM   
%f_max=Discrete_F(length(Discrete_F))*ones(1,M);% maximum frequencies for each VM
%PowerFunc=pow_p(F,3);
omega=1*ones(1,M);% cost portion of VM consumption
%P_Idle=0.5*ones(M,1);% Watt: Idle Power for each VM
E_max=100*ones(1,M);% Joule-maximum tolarated energy for each VM
%Delta=0.1*ones(1,M); %seconds-considered time of processing for each VM
W=1*ones(1,M);% MHZ frequency needed for communication channel
Zeta=[0.2;0.5]; %milli Watt
Zeta1=0.2; %milli Watt
Zeta2=0.5; %milli Watt
Zeta3=[0.5:0.25:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.25
Zeta4=[0.5:0.5:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.5

Omega=[0.2;0.5]; %milli Watt
Omega1=0.2; %milli Watt
Omega2=0.5; %milli Watt
Omega3=[0.5:0.25:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.25
Omega4=[0.5:0.5:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.5

%L_tot=60 + (80-60).*rand(1,1000); %carico uniformemente distribuito   70+-10 with avg(L_{tot})=70 (PMR=1.14285) %UTILIZZATO PER Scenario 2
%save('workload1000M.mat','L_tot');
% L_tot_Vector = load('workload1000M','L_tot');
% L_tot=L_tot_Vector.L_tot;
L_tot_Vector  = load('workload1000');
L_tot=L_tot_Vector.workload1000;
%L_tot = 8 + PMR .* rand(1,N);%task-size of each VM
f_opt_wokloads=zeros(N,M); % optimal frequencies of M VMs for each workload
%L_opt_wokloads=zeros(N,M); % optimal Load of M VMs for each workload
t_opt_wokloads=zeros(N,Q+1); % optimal rate of M VMs for each workload
%I=10;% number of window-slots
Opt_Comp_Energy=zeros(length(L_tot),M);
Opt_Comp_Energy1=zeros(length(L_tot),M);
Opt_reconfig_Energy=zeros(length(L_tot),M);
Opt_reconfig_Energy11=zeros(length(L_tot),M);
Opt_reconfig_Energy12=zeros(length(L_tot),M);
Opt_reconfig_Energy13=zeros(length(L_tot),M);
Opt_Comm_Energy=zeros(length(L_tot),M);
Opt_Comm_Energy1=zeros(length(L_tot),M);
Opt_Comm_Energy2=zeros(length(L_tot),M);
Opt_Comm_Energy3=zeros(length(L_tot),M);
Opt_Comm_Energy4=zeros(length(L_tot),M);

total_energy=zeros(length(L_tot),M);
total_energy1=zeros(length(L_tot),M);
total_energy2=zeros(length(L_tot),M);
total_energy3=zeros(length(L_tot),M);
total_energy4=zeros(length(L_tot),M);


%%%%%%% CVX solution
Opt_reconfig_Energy1=zeros(length(L_tot),1);
%MM=[20,30,40];
%for j=5:M,% Vm number
% total_energy1=zeros(length(T_vector),length(R_tot_vector),length(L_tot),M);
% total_energy2=zeros(length(T_vector),length(R_tot_vector),length(L_tot),M);
% total_energy3=zeros(length(T_vector),length(R_tot_vector),length(L_tot),M);
% total_energy4=zeros(length(T_vector),length(R_tot_vector),length(L_tot),M);

total_energy1=zeros(2,length(L_tot),M);
total_energy2=zeros(2,length(L_tot),M);


tic
%for tt=1:length(T_vector),
 %   T=T_vector(tt);
%for rr=1:length(R_tot_vector),
 %   R_tot=R_tot_vector(rr);
  %  R_tot_vector
%     for jj=1:3,
%     j=MM(jj);
%      R_tot1=R_tot.*ones(1,j);

 for jj=1:2
     if (jj==1)
     Discrete_F=Discrete_F2;
     else
     Discrete_F=Discrete_F1;
     end
 for j=1:M,
     
     if (jj==2) && (j<3)
             continue
         end
    %j=MM(jj);
    temp=0;
    f_zero=0;
    F=ones(j,1)*Discrete_F; % Discrete-frequency vector for M vms (Hemogenous VMs)
    PowerFunc=pow_p(F,3);
    P_Idle=0.5*ones(j,1);% Watt: Idle Power for each VM channel
    %f_max=0.933*ones(j,1);% maximum frequencies for each VM
    f_max=max(Discrete_F).*ones(1,j);
    Delta_F=zeros(j,1);
    for i=1:(length(L_tot)),
        
        %L_tot=L_tot(1+10*(i-1):1:10*i);
        %L_tot=1.867*5;
        %L_tot=200;
        %for t=1:I,% for window-slots for each workload
        %L_tot=10;
       % fprintf(1, 'L_tot value: %8.4e\n', L_tot(i));
        %fprintf(1, 'Number of L_tot \n');
        %disp(i);
        %fprintf('VM number= \n');
        %disp(j);
        %*******************Feasibility condition************************************

%         condizione_feas_workloadcomputation= (L_tot(i)<= sum(f_max(:).*T(1)));
%         condizione_feas_channel=(L_tot(i)<= R_tot*Q*(T_t-T(1))/2);
         
        condizione_feas_workloadcomputation= (L_tot(i)<= sum(f_max(:).*T(1)));
        condizione_feas_channel=(L_tot(i)<= R_tot*Q*(T_t-T(1))/2);
%        condizione_feas_channel=(L_tot(i)<= sum(R_tot1)*Q*(T_t-T)/2);

        
              
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
%    expression f_zero(I);% starting frequncy of each VM
    %f_zero=temp;
    %expression L_b(I);% background wokload for each VM before starting processing
    %expression temp_f(I);% temporary freuqncy for each VM while processing
   %expressions temp1(I) temp2(I);
   %expressions temp_ss(M,Q+1) t_opt_temp(M,Q) t_opt_sorted(M,Q+1);
    %L_b=zeros(1,M);
    %f_zero(:)=zeros(1,M);
    %variable R(j,I); % optimum rate of channels for each VM
    %variable f(I); % optimum frequency for each VM
    variable t_opt(j,Q+1); % optimum load for each VM
    %dual variables A B C D E
    %dual variable B(M)
    %t_opt_temp=t_opt(:, 2:Q+1);

    % objective optimization
    %minimize((omega(1)*E_max(1)*sum_square(vec(f)/f_max(1)))+k_e(2)*sum_square((vec(f)-vec(f_zero)))+(T_tot-Delta(1))*Zeta(2)*sum(power(2,2.*L(:)/(T_tot-Delta(1))*W(1))-1));
    %minimize (omega(1)*sum(pow_p(F,3)).*t_opt(:)+k_e(2)*sum((F-f_zero(:)).^2));
    %(omega(1)*sum(PowerFunc(:)'*t_opt(:)))  for Computaion: OK
    %(T_t-T(1))*Zeta(2)*sum(power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)%for Communication: OK
    %minimize ((omega(1)*sum(PowerFunc(:).*t_opt(:)))+((T_t-T(1))*((Omega(2)*sum(((2.*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)))));
    minimize ((omega(1)*sum(PowerFunc(:).*t_opt(:)))+((T_t-T(1))*((Zeta(2)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)))));
    
    %s.t. constrainsts
    subject to
    sum(sum(F(:,2:Q+1).*t_opt(:,2:Q+1)))-L_tot(i)==0;
    for  kk=1:j
    sum(t_opt(kk,:))-T(1)==0;
    end
    0<= t_opt(:) <= T(1);
    sum((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1)))-R_tot<=0;
    
    %B: sum(t_opt(1,:))-T(1)==0;
    %C: sum(t_opt(2,:))-T(1)==0;
    %D: 0<= F(:) <= f_max(1);
    %E: 0<= t_opt_temp <= T(1);
    %C: -L<=0;
    %F: -R<=0;
    
        %temp_f=f;
%     for k=1:I
%     temp1(k)=f(k);
%     temp2(k)=t_opt(k);
%     temp3(k)=temp1(k)+temp2(k);
%     end
%     for kk=1:I
%     sum(temp3)-L_tot(kk)==0;
%     end
    cvx_end
    %Reconfiguration cost******begin
%     for  mm=1:M
%         [t_opt_sorted(mm,:),t_opt_index]=sort(t_opt(mm,:),'descend');
%     end
%     test=issorted(t_opt_sorted, 'rows');
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
    Opt_reconfig_Energy1(i)=sum(Delta_F(:));
    Delta_F=zeros(j,1);
    %Reconfiguration cost******end
    %%%%%%%%%%%%%%%%%%%%%%%%Results
    Opt_Comp_Energy(i,j)=omega(1)*sum(PowerFunc(:).*t_opt(:));
    %Opt_reconfig_Energy1(i,j)=k_e(1)*sum(Opt_reconfig_Energy1(i));
     Opt_reconfig_Energy11(i,j)=k_e(1)*sum(Opt_reconfig_Energy1(i));
     Opt_reconfig_Energy12(i,j)=k_e(2)*sum(Opt_reconfig_Energy1(i));
%     Opt_reconfig_Energy13(i,j)=k_e(3)*sum(Opt_reconfig_Energy1(i));
    %Opt_Comm_Energy(i,j)=(T_t-T(1))*((Zeta(2)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
    % scenario 1
    Opt_Comm_Energy1(i,j)=(T_t-T(1))*((Zeta1*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
    Opt_Comm_Energy2(i,j)=(T_t-T(1))*((Zeta2*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
    Opt_Comm_Energy3(i,j)=(T_t-T(1))*((Zeta3(j)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
    Opt_Comm_Energy4(i,j)=(T_t-T(1))*((Zeta4(j)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
%     % scenario 2
    %Opt_Comm_Energy(i,j)=(T_t-T(1))*((Omega(2)*sum(((2*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)));
%     Opt_Comm_Energy1(i,j)=(T_t-T(1))*((Omega1*sum(((2*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)));
%     Opt_Comm_Energy2(i,j)=(T_t-T(1))*((Omega2*sum(((2*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)));
%     Opt_Comm_Energy3(i,j)=(T_t-T(1))*((Omega3(j)*sum(((2*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)));
%     Opt_Comm_Energy4(i,j)=(T_t-T(1))*((Omega4(j)*sum(((2*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)));
    %total_energy(i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy(i,j)+Opt_Comm_Energy(i,j);
    
%     total_energy1(i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy1(i,j);
%     total_energy2(i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy2(i,j);
%     total_energy3(i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy3(i,j);
%     total_energy4(i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy4(i,j);

     total_energy1(jj,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy11(i,j)+Opt_Comm_Energy2(i,j);
     total_energy2(jj,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy12(i,j)+Opt_Comm_Energy2(i,j);

 
%     total_energy1(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy1(i,j);
%     total_energy2(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy2(i,j);
%     total_energy3(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy3(i,j);
%     total_energy4(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy4(i,j);
%     
     %cvx_cputime_1(i,j)=toc;
    cvx_cputime_2(i,j)=cvx_cputime;
    %f_opt_wokloads(i,:)=f(1);
    %L_opt_wokloads(i,:)=L(1);
    t_opt_wokloads(i,:)=t_opt(j,:);
    
    
    
    
    
    
%     temp=f;
%     if (i==length(L_tot_Vector)/I)
%         temp=0;% starting frequncy of each VM in new loop
%     end
    % end % for time-slots
    end % for jobs (workloads)
 end % for VMs
 end % for various k_e
%end % for R_tot_vector
%end % for T_vector
time=toc;
%save('Myapproach_energy1000M2.mat','total_energy','Opt_Comp_Energy','Opt_reconfig_Energy','Opt_Comm_Energy','cvx_cputime_2');    
%%%%%%% simulation plots

figure(100001)
plot([1:1:M],sum(total_energy1)/length(L_tot),'--+',[1:1:M],sum(total_energy2)/length(L_tot),'--*',[1:1:M],sum(total_energy3)/length(L_tot),'--o',[1:1:M],sum(total_energy4)/length(L_tot),'--+');
xlabel('VM');
ylabel('Average E_{total} (Joule)');
legend('\zeta=0.2 (mW)','\zeta=0.5 (mW)','\zeta=0.5+0.25(i-1) (mW)','\zeta=0.5+0.5(i-1) (mW)');
grid on

% figure(100002)
% plot([1:1:M],sum(total_energy1)/length(L_tot),'--+',[1:1:M],sum(total_energy2)/length(L_tot),'--*',[1:1:M],sum(total_energy3)/length(L_tot),'--o');
% xlabel('VM');
% ylabel('Average E_{CPU} (Joule)');
% legend('k_e=0.005 (Joule/(MHz)^2)','k_e=0.05 (Joule/(MHz)^2)','k_e=0.1 (Joule/(MHz)^2)');
% grid on
% 
% figure(100003)
% plot([1:1:M],sum(total_energy1)/length(L_tot),'--+',[1:1:M],sum(total_energy2)/length(L_tot),'--*',[1:1:M],sum(total_energy3)/length(L_tot),'--o');
% xlabel('VM');
% ylabel('Average E_{Reconf} (Joule)');
% legend('k_e=0.005 (Joule/(MHz)^2)','k_e=0.05 (Joule/(MHz)^2)','k_e=0.1 (Joule/(MHz)^2)');
% grid on
% 
% figure(100004)
% plot([1:1:M],sum(total_energy1)/length(L_tot),'--+',[1:1:M],sum(total_energy2)/length(L_tot),'--*',[1:1:M],sum(total_energy3)/length(L_tot),'--o');
% xlabel('VM');
% ylabel('Average E_{Net} (Joule)');
% legend('k_e=0.005 (Joule/(MHz)^2)','k_e=0.05 (Joule/(MHz)^2)','k_e=0.1 (Joule/(MHz)^2)');
% grid on

% figure(1)
% plot([1:1:1000],cvx_cputime_1(:,2),'--+',[1:1:1000],cvx_cputime_1(:,10),'--*');
% xlabel('Workload');
% ylabel('Run-time(seconds)');
% legend('VM=2','VM=10');
% grid on

% test=sum(t_opt_wokloads)/1000;
% figure(6)
% %plot([1:1:100],sum(test(:)),'--+');
% xlabel('Discrete Frequencies');
% ylabel('Average consumed time for each ADF');
% legend('VM=10');
% grid on



figure(1)
plot([1:1:N],f_opt_wokloads(:,2),'--+',[1:1:N],f_opt_wokloads(:,M),'--*');
xlabel('Workloads');
ylabel('optimum frequency for each VM');
legend('f-opt for VM=2','f-opt for VM=Last');
grid on

figure(2)
plot([1:1:M],f_opt_wokloads(1,:),'--+',[1:1:M],f_opt_wokloads(N,:),'--*');
xlabel('VMs');
ylabel('optimum frequency for each VM');
legend('f-opt:first workload','f-opt: last workload');
grid on

figure(3)
plot([1:1:N],R_opt_wokloads(:,2),'--+',[1:1:N],R_opt_wokloads(:,M),'--*');
xlabel('Workloads');
ylabel('optimum channel rate for each VM');
legend('R\_* for VM=2','R\_* for VM=Last');
grid on

figure(4)
plot([1:1:N],f_opt_wokloads(:,2),'--+',[1:1:N],f_opt_wokloads(:,M),'--*');
xlabel('Workloads');
ylabel('optimum frequency for each Workload');
legend('f-opt for VM=2','f-opt for VM=Last');
grid on

figure(5)
plot([1:1:N],Opt_Comp_Energy(:,M),'--+',[1:1:N],Opt_reconfig_Energy(:,M),'--o',[1:1:N],Opt_Comm_Energy(:,M),'--*')
xlabel('Workloads');
ylabel('optimum computation/reconfiguration/communication costs for each Workload');
legend('Computation Cost for M=last','reconfiguration Cost for M=last','Communication Cost for M=last');
grid on

figure(6)
plot([1:1:M],sum(total_energy(:,M))/N,'--+')
xlabel('VMs');
ylabel('Average total consumed energy for each VMs');
legend('Avarage total energy consumed');
grid on

