clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%difinite variables
M=10000;% numbrs of virtual machines
M_vector=[5000:1000:10000];
N=10;% workload numbers
Q=3;% numbers of Discrete ranges for Frequency of each VM
%T=0.1*ones(1,M);% seconds/ time for each Vm for computation
T=3;
T_t=5;
%T_vector=[1:0.1:T_t-1];
R_tot=1000000;% Mbit / seconds
%R_tot_vector=[10,100];
%PMR=1.25;
temp=0;
ALPHA=[0.1,0.5,0.9];

k_e=[0.005]*100;% Joule/(Mega Hz)^2 % scenario 1
%Discrete_F=[0.15, 1.867, 2.133,2.533,2.668];%(GhZ)% for Interl Nehalem Quad-core Processor in INFOCOM 2009 G.V. Laszewski et al
%Discrete_F=[0.3, 0.533, 0.667, 0.800, 0.933];%(GhZ)% for power-scalable real cluster Crusoe with CPU type TM-5800 in CLUSTER 2006 in http://dx.doi.org/10.1109/CLUSTR.2006.311839
%F=ones(M,1)*Discrete_F; % Discrete-frequency vector for M vms (Hemogenous VMs)
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


%F1=ones(M,1)*Discrete_F(2:length(Discrete_F)); % Discrete-frequency vector for M vms (Hemogenous VMs) without F_0 which is belongs to Idle state
%Delta_F=zeros(M,1);  % (f(new)-f(old))^2
%f_zero=zeros(M,1);% starting frequncy of each VM   
%f_max=Discrete_F(length(Discrete_F))*ones(1,M);% maximum frequencies for each VM
%PowerFunc=pow_p(F,3);
omega=1;% cost portion of VM consumption
%P_Idle=0.5*ones(M,1);% Watt: Idle Power for each VM
%E_max=100*ones(1,M);% Joule-maximum tolarated energy for each VM
%Delta=0.1*ones(1,M); %seconds-considered time of processing for each VM
W=1;% MHZ frequency needed for communication channel

Omega1=5*10^(-3);%Watt
Omega2=50*10^(-3);%Watt

% Omega=[0.2;0.5]; %milli Watt
% Omega1=0.2; %milli Watt
% Omega2=0.5; %milli Watt
% Omega3=[0.5:0.25:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.25
% Omega4=[0.5:0.5:100]; %milli Watt for 10 VM (increament and depend on VM) with slip 0.5

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
%mL_tot=900;
mL_tot=90000;
pmr=1.5;
a=mL_tot*(pmr-1);
if pmr>1
%L_tot=(mL_tot-a) + 2.*(a).*rand(1,N);
%save('Ltot900PMR1.5.mat','L_tot');
%save('Ltot90000PMR1.5.mat','L_tot');
%L_tot=load('Ltot900PMR1.5.mat','L_tot');
L_tot=load('Ltot90000PMR1.5.mat','L_tot');
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
%I=10;% number of window-slots
% Opt_Comp_Energy=zeros(length(L_tot),length);
% Opt_Comp_Energy1=zeros(length(L_tot),M);
% Opt_reconfig_Energy=zeros(length(L_tot),M);
% Opt_reconfig_Energy11=zeros(length(L_tot),M);
% Opt_reconfig_Energy12=zeros(length(L_tot),M);
% Opt_reconfig_Energy13=zeros(length(L_tot),M);
% Opt_Comm_Energy=zeros(length(L_tot),M);
% Opt_Comm_Energy1=zeros(length(L_tot),M);
% Opt_Comm_Energy2=zeros(length(L_tot),M);
% Opt_Comm_Energy3=zeros(length(L_tot),M);
% Opt_Comm_Energy4=zeros(length(L_tot),M);
% 
% total_energy=zeros(length(L_tot),M);
% total_energy1=zeros(length(L_tot),M);
% total_energy2=zeros(length(L_tot),M);
% total_energy3=zeros(length(L_tot),M);
% total_energy4=zeros(length(L_tot),M);


%%%%%%% CVX solution
%Opt_reconfig_Energy1=zeros(length(L_tot),1);
MM=M_vector;
Opt_Comp_Energy=zeros(length(L_tot),length(MM));
Opt_Comp_Energy1=zeros(length(L_tot),length(MM));

Opt_reconfig_Energy=zeros(length(L_tot),length(MM));
Opt_reconfig_Energy11=zeros(length(L_tot),length(MM));

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
keep_NaNslot=zeros(length(L_tot),length(MM));
power_calc=zeros(length(L_tot),length(MM));
tic
%for tt=1:length(T_vector),
%   T=T_vector(tt);
% for rr=1:length(R_tot_vector),
%     R_tot=R_tot_vector(rr);
%     R_tot_vector
%     for jj=1:3,
%     j=MM(jj);
  %   R_tot1=R_tot.*ones(1,j);
 Total_Time=zeros(length(L_tot),M,Q+1);
 t_opt_wokloads_per_Q=zeros(1,Q+1);
 %%
    %%code body
 %for j=5:M,
     Opt_reconfig_Energy1=zeros(length(L_tot),length(MM));
Omega3_vector=[45:-5:10].*10^(-3);
  for jj=1:length(MM),
    j=MM(jj);
    Total_Time=zeros(length(L_tot),j,Q+1);
    temp=0;
    f_zero=0;
    F=ones(j,1)*Discrete_F; % Discrete-frequency vector for M vms (Hemogenous VMs)
    PowerFunc=pow_p(F,3);
    P_Idle=0.0005*ones(j,1);% Watt: Idle Power for each VM channel
    %f_max=0.933*ones(j,1);% maximum frequencies for each VM
    f_max=max(Discrete_F).*ones(1,j);
    Delta_F=zeros(j,1);
    
    Omega3=(5+45*rand(1,1)).*10^(-3);% milli watt random variable
    Omega3=Omega3_vector(jj);% milli watt random variable

    
    for i=1:(length(L_tot)),
        tic
       if L_tot(i)> sum(f_max(:).*T(1))
          error('Problem SLA violation')
       end
       
       if L_tot(i)==0
           
          continue
       end
                       
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

         condizione_feas_workloadcomputation= (L_tot(i)<= sum(f_max(:).*T(1)));
         condizione_feas_channel=(L_tot(i)<= R_tot*Q*(T_t-T(1)/2));
         
%        condizione_feas_workloadcomputation= (L_tot(i)<= sum(f_max(:).*T));
%        condizione_feas_channel=(L_tot(i)<= R_tot1*Q*(T_t-T)/2);
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
    variable t_opt(j,Q+1) nonnegative; % optimum load for each VM
    %dual variables A B C D E
    %dual variable B(M)
    %t_opt_temp=t_opt(:, 2:Q+1);

    % objective optimization
    %minimize((omega(1)*E_max(1)*sum_square(vec(f)/f_max(1)))+k_e(2)*sum_square((vec(f)-vec(f_zero)))+(T_tot-Delta(1))*Zeta(2)*sum(power(2,2.*L(:)/(T_tot-Delta(1))*W(1))-1));
    %minimize (omega(1)*sum(pow_p(F,3)).*t_opt(:)+k_e(2)*sum((F-f_zero(:)).^2));
    %(omega(1)*sum(PowerFunc(:)'*t_opt(:)))  for Computaion: OK
    %(T_t-T(1))*Zeta(2)*sum(power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)%for Communication: OK
    %minimize ((omega(1)*sum(PowerFunc(:).*t_opt(:)))+((T_t-T(1))*((Omega(2)*sum(((2.*(sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/((T_t-T(1))*W(1))).^(1/ALPHA(2)))))+sum(P_Idle(:)))));
    minimize ((omega(1)*sum(PowerFunc(:).*t_opt(:)))+((T_t-T(1))*((Omega1*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)))));
    
    %s.t. constrainsts
    subject to
    sum(sum(F(:,2:Q+1).*t_opt(:,2:Q+1)))-L_tot(i)==0;
    for  kk=1:j
    0<=sum(t_opt(kk,:))-T(1)<=0;
   % (2*(T_t-T(1))).*(sum(F(kk,:).*t_opt(kk,:))./(sum(F(kk,:).*t_opt(kk,:),2)))+T(1)-T_t<=0;
    end
    %for  kk=1:j
     %   for  fff=1:Q+1
            0<= t_opt <= T(1);
      %  end
    %end
    0<= F <= max(Discrete_F);

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
    Opt_reconfig_Energy1(i,jj)=sum(Delta_F(:));
    Delta_F=zeros(j,1);
    %%
    %%Results
    for titr=1:j
        for titr1=1:Q+1
            if ~isnan(t_opt(titr,titr1))
                power_calc(titr,titr1)=PowerFunc(titr,titr1)*t_opt(titr,titr1);
            else
                power_calc(titr,titr1)=0;
            end
        end
    end
     Opt_Comp_Energy(i,jj)=omega(1)*sum(power_calc(:));
    % CPU_energy1(tt,i,j)=Opt_Comp_Energy(i,j);
      Opt_reconfig_Energy11(i,jj)=k_e(1)*Opt_reconfig_Energy1(i,jj);
     
      
      %Opt_Comm_Energy(i,j)=(T_t-T(1))*((Zeta(2)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T(1))*W(1)))-1)))+sum(P_Idle(:)));
    RTT=70;
    % scenario 1
    Opt_Comm_Energy1(i,jj)=(sum(Omega1*(RTT.*2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T))).^2)+sum(P_Idle(:)));
    Opt_Comm_Energy2(i,jj)=(sum(Omega2*(RTT.*2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T))).^2)+sum(P_Idle(:)));
    Opt_Comm_Energy3(i,jj)=(sum(Omega3*(RTT.*2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T))).^2)+sum(P_Idle(:)));
    %Opt_Comm_Energy1(i,jj)=(T_t-T(1))*((Zeta5*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
    %Opt_Comm_Energy2(i,jj)=(T_t-T(1))*((Zeta5*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
    %Opt_Comm_Energy3(i,jj)=(T_t-T(1))*((Zeta5*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
    
    total_energy1(i,jj)=Opt_Comp_Energy(i,jj)+Opt_reconfig_Energy11(i,jj)+Opt_Comm_Energy1(i,jj);
    total_energy2(i,jj)=Opt_Comp_Energy(i,jj)+Opt_reconfig_Energy11(i,jj)+Opt_Comm_Energy2(i,jj);
    total_energy3(i,jj)=Opt_Comp_Energy(i,jj)+Opt_reconfig_Energy11(i,jj)+Opt_Comm_Energy3(i,jj);

    % scenario 11
%     Opt_Comm_Energy1(i,j)=(T_t-T(1))*((Zeta1*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
%     Opt_Comm_Energy2(i,j)=(T_t-T(1))*((Zeta2*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
%     Opt_Comm_Energy3(i,j)=(T_t-T(1))*((Zeta3(j)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
%     Opt_Comm_Energy4(i,j)=(T_t-T(1))*((Zeta4(j)*sum((power(2,2.*((sum(F(:,2:Q+1).*t_opt(:,2:Q+1),2))/(T_t-T)*W(1)))-1)))+sum(P_Idle(:)));
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
%     total_energy1(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy1(i,j);
%     total_energy2(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy2(i,j);
%     total_energy3(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy3(i,j);
%     total_energy4(tt,rr,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy4(i,j);
    
 %   total_energy1(tt,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy1(i,j);
 %   total_energy2(tt,i,j)=Opt_Comp_Energy1(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy2(i,j);
    %total_energy3(tt,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy3(i,j);
    %total_energy4(tt,i,j)=Opt_Comp_Energy(i,j)+Opt_reconfig_Energy1(i,j)+Opt_Comm_Energy4(i,j);

    
     %cvx_cputime_1(i,j)=toc;
    cvx_cputime_2(i,j)=cvx_cputime;
    %f_opt_wokloads(i,:)=f(1);
    %L_opt_wokloads(i,:)=L(1);
   % if (t_opt~=NaN)
    t_opt_wokloads(i,:)=t_opt(j,:);
    
    if isnan(t_opt)
    keep_NaNslot(i,jj)=1;
    end
    t_opt_wokloads(any(isnan(t_opt_wokloads),2),:) = [];
   % t_opt_wokloads(~any(isnan(t_opt_wokloads)));
   % t_opt_wokloads(isnan(t_opt_wokloads))=0;
   %end
    Total_Time(i,:,:)=t_opt;
    

    end % for jobs (workloads)
  end % for VMs
  
  
%end % for R_tot_vector
%end % for T_vector
time=toc;

 t_opt_wokloads(any(isnan(t_opt_wokloads),2),:) = [];
% for j=1:Q+1
% t_opt_wokloads_per_Q(j)=mean(t_opt_wokloads(:,j));
% end
% bar(t_opt_wokloads_per_Q/sum(t_opt_wokloads_per_Q));
% title(['mean(L_{tot})=' num2str(mean(L_tot)) ', UNI2\_IMC2010\_DC2, max(L_{tot})= ' num2str(max(L_tot)) ', M= ' num2str(M)...
%     ', T= ' num2str(T) ', N= ' num2str(N) ', PMR= ' num2str(pmr) ', F_Q= ' num2str(max(Discrete_F))]);
% grid on
total_energyt1=zeros(1,length(MM));
total_energyt2=zeros(1,length(MM));
total_energyt3=zeros(1,length(MM));

Opt_Comp_Energyt1=zeros(1,length(MM));

Opt_Comm_Energyt1=zeros(1,length(MM));
Opt_Comm_Energyt2=zeros(1,length(MM));
Opt_Comm_Energyt3=zeros(1,length(MM));

tot_reconfig_Energy1=zeros(1,length(MM));



for jj=1:length(MM)
    Opt_Comp_Energyt1(jj)=mean(Opt_Comp_Energy(:,jj));
    
    tot_reconfig_Energy1(jj)=mean(Opt_reconfig_Energy11(:,jj));
    
    Opt_Comm_Energyt1(jj)=mean(Opt_Comm_Energy1(:,jj));
    Opt_Comm_Energyt2(jj)=mean(Opt_Comm_Energy2(:,jj));
    Opt_Comm_Energyt3(jj)=mean(Opt_Comm_Energy3(:,jj));
  

    total_energyt1(jj)=mean(total_energy1(:,jj));
    total_energyt2(jj)=mean(total_energy2(:,jj));
    total_energyt3(jj)=mean(total_energy3(:,jj));

end

save('Myapproach_energyScnerioN100F1Big.mat','total_energyt1','total_energyt2','total_energyt3','tot_reconfig_Energy1',...
     'Opt_Comm_Energyt1','Opt_Comm_Energyt2','Opt_Comp_Energyt1','Opt_Comm_Energyt3','cvx_cputime_2'); 
%%
%%%%%%% simulation plots
figure(2000) %fig 4a
plot(1:1:length(MM),total_energyt1,'-b',1:1:length(MM),total_energyt2,'-k',1:1:length(MM),total_energyt3,'-r');
xlabel('$VMs$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}_{tot}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\Omega_i=5\:(mW)$','$\Omega_i=50\:(mW)$','$\Omega_i=[5,50]\:(mW)$','Interpreter','latex','FontSize',40);
grid on

figure(2001) %fig 4b
plot(1:1:length(MM),Opt_Comm_Energyt1,'-b',1:1:length(MM),Opt_Comm_Energyt2,'-k',1:1:length(MM),Opt_Comm_Energyt3,'-r');
xlabel('$VMs$','Interpreter','latex','FontSize',40);
ylabel('$\overline{\mathcal{E}}^{CMc}\:(Joule)$','Interpreter','latex','FontSize',40);
legend('$\Omega_i=5\:(mW)$','$\Omega_i=50\:(mW)$','$\Omega_i=[5,50]\:(mW)$','Interpreter','latex','FontSize',40);
grid on
