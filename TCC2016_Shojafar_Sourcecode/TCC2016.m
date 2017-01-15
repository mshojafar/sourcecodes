
clear all;
close all;
clc;
%% Loading and assignments

load('X');
load('EPS');
load('DELTA_IP');

%% %PARAMETRI GENERALI
n=2000;                                                                    %numero of slot
nclient=1;                                                                 % number of clients
Ts = 2;                                                                    %time slot (sec)
kappa=10000;
up_th=0.70;
% X=rand(1,2*n);
% save('X');
% EPS=rand(1,2*n);
% save('EPS');
%DELTA_IP=rand(1,n);
%save('DELTA_IP');


%Ew =(-1.732+25)+(2*1.732*EW(1:N));
%Ew =(7)+(EW(1:n));
E_W=zeros(nclient,n);                                                       %TCP/IP energy consumption for each slot
%E_TOT=zeros(nclient,n);                                                         %TOT energy consumption for each slot
Mu=zeros(nclient,n);  %[ones(nclient,1) zeros(nclient,n)];                                      %moltiplicatore
E_TOT=zeros(1,n);                                                         %TOT energy consumption for each slot


avg_r=zeros(nclient,n);
q=zeros(nclient,n);                                                        %service rate OUT-queue sistema
s=zeros(nclient,n);                                                        %service rate IN-queue sistema
a=zeros(nclient,n);                                                        %incoming rate IN-queue sistema
Jobs=a;
r=zeros(nclient,n);                                                        %arrival rate in INPUT queue
L_R=zeros(nclient,n);                                                      %real workload incoming to the buffer/queue
%L_tot_vect = 10+ (200) * rand(1,n);
%% PMR=1.25, L_tot=8+-2   L_tot=12+-3
%save('workloads8PMR125No100000Journal.mat', 'L_tot_vect');                  %loading 100 fixed random workloads
L_tot_vect=load('workloads8PMR125No100000Journal.mat','L_tot_vect');                %loading 2000 fixed random workloads
L_tot_vect = L_tot_vect.L_tot_vect;
%save('workloads100PMR2No20Journal.mat','L_tot_vect'); %%loading 1000 fixed random workloads-test
%L_tot_vect=load('workloads50PMR1.5No1000Journal.mat','L_tot_vect'); %%loading 1000 fixed random workloads-test
%L_tot_vect=load('workloads100PMR2No20Journal.mat','L_tot_vect'); %%loading 1000 fixed random workloads-test
%L_tot_vect = L_tot_vect.L_tot_vect;
%L_tot_vect = load('workload1M_mds1_MSR.mat','L_tot'); % real MSR workload
%L_tot_vect = load('per_second_request_day_50.mat','per_second_request'); % real WorldCUp 98 50th day workload: 24*60*60 in seconds

%L_tot_vect = L_tot_vect.per_second_request;
L_tot_vect=L_tot_vect(1:n);
% L_tot_vect(2)=14;
% L_tot_vect(3)=194;
% L_tot_vect(4)=10;
% L_tot_vect(5)=40;
% L_tot_vect(11)=4.5;
% L_tot_vect(20)=14;
% L_tot_vect(30)=150;
% L_tot_vect(40)=10;
% L_tot_vect(45)=400;
% L_tot_vect(48)=1500;
%L_tot_vect(300)=190;
%L_tot_vect(400)=3;
%%
%
%   for x = 1:10
%       disp(x)
%   end
%
%L_tot_vect(100)=10;
%%
%
% $$e^{\pi i} + 1 = 0$$
%

a=L_tot_vect;

L_tot=zeros(1,n);                                                          %workload per slot of the sistema
lambda=zeros(nclient,n);                                                   %arrival rate in INPUT queue

%% Critical parameter%s
energy_ave=0.125;
%energy_ave_vector=[20, 10,1,0.5,0.25,0.125,0.0625,0.03125,0.01, 0.005];   % first test
%energy_ave_vector=[10, 0.25];   % second test
%energy_ave_vector=[0.75];   % third test
N_O=240;       %old=240000                                                         %MAXIMUM CAPACITY of OUT-queue(byte)
N_I=240;       %old=240000                                                           %MAXIMUM CAPACITY of IN-queue(byte)
%r_min=3.5; %we suppose average of a=8                                                             %(byte/slot)
r_max=16;                                                                %(byte/slot)
%a_min=r_min;
%a_max=r_max;

%% PARAMETRI DI CANALE
DELTA_IP_MAX = 10;                                                         %max ip layer delay = 10*Ts(time slot)
%DELTA_IP1 = (round(DELTA_IP_MAX * DELTA_IP(1:n)));
DELTA_IP1 = ((DELTA_IP_MAX * DELTA_IP(1:n)));

RTT = zeros(nclient,n);
b = 2;
MSS = 120;
beta=1;
%teta=0.5;
sigma = zeros(nclient,n);
v = 13;                                                                    %speed m/s (13 m/s is the v of the paper: h=0.95)
h = 0.82^(v*Ts/100);                                                       %correlation coefficient

A = 90.2514;
B = 3.4998;
C = 10^0.10942;
% A = 274.7229;
% B = 7.9932;
% C = 10^(-0.15331);
x = X(1:2*n);
esp = EPS(1:2*n);
for k = 2:2*n
    x(1,k) = (h^0.5)*x(1,k-1) + ((1-h)^0.5)*esp(1,k);
end
x = x(n+1:2*n);
% x = x - mean(x);
x = x - min(x);
% x = x - mean(x);
x = (1/max(x)) * x;%normalization
x = -3^0.5 + 2*(3^0.5)*x;


a_0 = exp(-0.5*(0.1*log(10))^2);
z = a_0*10.^(0.1*x);
%Pl = (C + (A/((C*B)^2)) * gammainc(1,C*B)) * (z./E_W);
% Pl = (0.7*10^-2)*ones(1,N);


K_0 = 1*( ((3/(2*b))^0.5) * MSS ) / (C+(A/(C*B^2)) * gammainc(1,C*B))^-0.5;

%r_min_vector=zeros(1,length(energy_ave_vector));
sigma_min_anal=K_0*((a_0*(10^(-0.1*sqrt(3))))^(1/2))/(DELTA_IP_MAX);
%r_min_vector=sigma_min.*sqrt(energy_ave_vector);
r_min_anal=sigma_min_anal.*sqrt(energy_ave);
%r_min_vector=[0.01, 0.1, 1, r_min_vector_anal, 4, 5];
% K_0 = 39124;

% K_0 = ( ((3/(2*b))^0.5) * MSS ) / (C+(A/(C*B^2)) * gamma(C*B))^-0.5;
% K_0 = ( ((3/(2*b))^0.5) * MSS ) / (C+((A * gammainc(1,C*B)^-0.5))/(C*B^2))
%
%% DVFS Block definitions

VM=[1:1:40];%numero Macchine Virtuali
NumSimulazioni=length(VM);

costo_tot_VM=zeros(1,NumSimulazioni);
costo_tot_VM_f_succ=zeros(1,NumSimulazioni);
costo_tot_VM_f_prec=zeros(1,NumSimulazioni);
costo_tot_VM_TS=zeros(1,NumSimulazioni);
costo_VM_TS_switch=zeros(1,NumSimulazioni);
costo_tot_VM_TSpiuSwitch=zeros(1,NumSimulazioni);


f_VM_prima=zeros(1,NumSimulazioni);
f_VM_ultima=zeros(1,NumSimulazioni);


iter_VM=zeros(1,NumSimulazioni);

L_VM_prima=zeros(1,NumSimulazioni);
L_VM_ultima=zeros(1,NumSimulazioni);

%SET UTILIZZATO PER FIG.2
%chi_completo=[0.5:0.25:100]; %1000 - BLU
%chi_completo=0.5.*ones(1,100);%2000 - VERDE
%chi_completo=0.2*ones(1,100);%3000 - ROSSO
%chi_completo=[0.5:0.5:100];%4000 - NERO

%chi_completo=[0.5:0.25:100];% heterogenous
chi_completo=0.5.*ones(1,1000);% hemogenous

%inizializzazione variabili

k_e=0.5; %0.005; %0.05; ke grande e piccolo
R_tot=200; %C_max=15; %(Mb/s)
T_tot=5; %(s)


%rand('twister',0);
%Jobs=8.*ones(1,1000);  %carico deterministico
%Jobs=6 + (10-6).*rand(1,1000); %carico uniformemente distribuito   8+-2   %UTILIZZATO PER FIG.2, FIG.6 e FIG.7
%Jobs=2 + (14-2).*rand(1,1000); %carico uniformemente distribuito   8+-6
%Jobs=(16).*rand(1,1000); %carico uniformemente distribuito   8+-8
%% Sigma/RTT calculation block

%sigma_vector(nclient,1)=[0.5, 1, 10, 100];
for i=1:n
    if (i==1)
        RTT(nclient,1)=0;
        sigma(nclient,1)=0;
        %sigma(nclient,2)=1;
    else
        RTT(nclient,i) = 0.75 * RTT(nclient,i-1) + 0.25 * DELTA_IP1(nclient,i);
        sigma(nclient,i) = ( K_0 * ((z(nclient,i))^0.5) ) / RTT(nclient,i);
    end
end
sigma_min_simulation=min(sigma(nclient,2:n));
RTT_max_simulation=max(RTT(nclient,:));

%% MINIBOOK MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test 5 parameters
q0_vector1=[0.01, 0.1, 1, 10, 20, 50, 100, 150, 200, N_O];  %third test
%q0_vector1=[0.01];  %second test
%q0_vector1=[240];  %third test
s0_vector1=[0.01];  %third test
%s0_vector1=[10^-2, N_I];  %third test
%teta_vector1=[0.01, 0.5, ]; %third test
teta_vector1=[0.5]; %8th test

J_P_vector=[0.25, 0.75];


% E_W_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% E_TOT_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% Muinf_vector=zeros(NumSimulazioni,length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% s_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% q_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% r_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% L_tot_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% lambda_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
% Utility_vector=zeros(NumSimulazioni, length(energy_ave_vector), length(q0_vector1), length(s0_vector1), length(teta_vector1));
%
%

E_W_vector=zeros(length(q0_vector1),length(J_P_vector));
E_TOT_vector=zeros(length(q0_vector1),length(J_P_vector));
Muinf_vector=zeros(length(q0_vector1),length(J_P_vector));
s_vector=zeros(length(q0_vector1),length(J_P_vector));
q_vector=zeros(length(q0_vector1),length(J_P_vector));
r_vector=zeros(length(q0_vector1),length(J_P_vector));
L_tot_vector=zeros(length(q0_vector1),length(J_P_vector));
lambda_vector=zeros(length(q0_vector1),length(J_P_vector));
Utility_vector=zeros(length(q0_vector1),length(J_P_vector));


% E_W_vector=zeros(1,length(q0_vector1));
% E_TOT_vector=zeros(1,length(q0_vector1));
% Muinf_vector=zeros(1,length(q0_vector1));
% s_vector=zeros(1,length(q0_vector1));
% q_vector=zeros(1,length(q0_vector1));
% r_vector=zeros(1,length(q0_vector1));
% L_tot_vector=zeros(1,length(q0_vector1));
% lambda_vector=zeros(1,length(q0_vector1));
% Utility_vector=zeros(1,length(q0_vector1));

% E_W_vector=zeros(1,length(r_min_vector));
% E_TOT_vector=zeros(1,length(r_min_vector));
% Muinf_vector=zeros(1,length(r_min_vector));
% s_vector=zeros(1,length(r_min_vector));
% q_vector=zeros(1,length(r_min_vector));
% r_vector=zeros(1,length(r_min_vector));
% L_tot_vector=zeros(1,length(r_min_vector));
% lambda_vector=zeros(1,length(r_min_vector));
% Utility_vector=zeros(1,length(r_min_vector));
%
% E_TOT_vector=zeros(1,NumSimulazioni);
W_range=zeros(1,n);
W_range_old=zeros(1,n);
W_range_estimate=zeros(1,n); % estimated to reach the up-limit
needed_new_VM=zeros(1,n);
mu_opt1=0;
SLA_violate=zeros(1,n);

% for j=1:M
%     if (Active_VM_rate(j)==1)
%         on_servers(1,j)=j;
%         l=l+1;
%         on_servers_list(j)=j;
%     end
% end

M1_vector=zeros(1,n);

for k=NumSimulazioni:NumSimulazioni
    tic
    k
    M=VM(k);
    %M1=min(M,ceil(mean(L_tot_vect)/10)+ceil((max(L_tot_vect)-mean(L_tot_vect))/10));
    %M1=min(M,ceil(mean(L_tot_vect)/10));
    M1=4;
    M1_start=M1;
    M1_vector(i)=M1;
    %M1=M/2; % maximum available VM can work
    %       Active_VM_rate = rand(1,M);                                                             % applications variable indicator
    Active_VM_rate = zeros(1,M);                                                             % we try to run half+1 of |VM in one system at first
    VM_inactive_vector= zeros(1,M);
    %     if floor(M/4)<M1 % for M>30 %pp=1:floor(M/2) % for low M<30
    %         pp_max=M1;
    %     else
    %         pp_max=floor(M/4);
    %     end
    %
    pp_max=M1;
    for pp=1:pp_max
        
        Active_VM_rate(pp)=1;
    end
    
    %      save('Active_VM_rate_value.mat','Active_VM_rate');
    % Active_VM_rate=load('Active_VM_rate_value.mat','Active_VM_rate');
    % Active_VM_rate=Active_VM_rate.Active_VM_rate;
    %Active_VM_rate(:)=(Active_VM_rate(:) >= 0.5);
    
    
    
    %INIZIALIZZAZIONE PARAMETRI %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ni_opt=zeros(1,M);
    f_opt=zeros(1,M);
    L_opt=zeros(1,M);
    mu_opt=0;
    iter_VM=zeros(1,M);
    alpha_step=0.01; %<==  N.B. Parametro utilizzato all'interno della funzione di tracking - settare in modo congruo
    V_step=0;
    
    f_max=10.*ones(1,M);  %(Mb/s)
    
    
    
    %f_zero=f_max;
    f_zero=zeros(1,M);  % (Mb/s)
    chi=chi_completo(1:M);
    W=ones(1,M);
    %N0=ones(1,M);
    %g=ones(1,M);
    Th=2.*log(2).*(chi./W);
    
    N_o_i=1;
    
    g=10^(3);
    sigma1=0.5*10^(-6);
    f_opt_ON=10^(-6)*f_max;
    
    
    E_max=40.*ones(1,M);  % (mJ)
    E_idle=1.51.*E_max;
    
    alpha=1;
    
    omega=1.*ones(1,M);
    Delta=1; %caso OMOGENEO: delta uguale per tutti i canali
    T_ON= 10^(-3).*Delta;
    
    L_max=f_max.*Delta;
    L_b=zeros(1,M);  % NB. il codice funziona bene per L_b=0, altrimenti aggiungere modifica su f_opt.
    
    temp=2.*k_e+(2.*E_max.*omega./(f_max.^2));
    alpha_zero=2.*k_e./temp;
    alpha_mu=Delta./temp;
    f_ibernazione=alpha_zero.*f_zero;
    
    
    %FREQUENZE DISCRETE INIZIALIZZAZIONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q=6; %numero frequenze discrete;
    f_discrete=zeros(M,Q);
    for conta=1:M
        f_discrete(conta,:)=[0:f_max(conta)./(Q-1):f_max(conta)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    test 2 parametrs
    %q0_vector1=[1, 5, 8, 10, 25, 50, 100];  %first test
    % q0_vector1=[1, 5, 8, 10, 25, 50, 100, 150, 200, 240];  %second test
    % E_W_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % Muinf_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % s_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % q_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % r_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % L_tot_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    % lambda_vector=zeros(length(energy_ave_vector), length(q0_vector1));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for test1=1:length(sigma_vector)
    %     test1
    %
    
    
    %r_min=r_min_vector(test10);
    %    for test1=1:length(energy_ave_vector)
    %        test1
    %
    r_min=r_min_anal;
    for test2=1:length(q0_vector1)
        q(nclient,1)=q0_vector1(test2);
        %q(nclient,1)=240;
        
        for test11=1:length(J_P_vector)
            J_P=J_P_vector(test11);
            %J_P=0.25;                                                       % jacoupson coefficeint
            
            
            %         if test2<=8
            %             kappa=10^4;
            %         else
            %             kappa=10^3;
            %         end
            %q(nclient,1)=100;
            %for test3=1:length(s0_vector1)
            %s(nclient,1)=s0_vector1(test3);
            s(nclient,1)=240;
            %s(nclient,1)=r_min;
            %      for test31=1:length(teta_vector1)
            %  teta=teta_vector1(test31);
            teta=0.5;
            
            avg_r_old=zeros(nclient,n);
            avg_r(nclient,1)=r_min;
            avg_r_old(nclient,1)=r_min;
            active_VM=zeros(1,n);                                             % numbers of active servers in each time slot
            active_VM_old=zeros(1,n);                                             % numbers of active servers in each time slot
            active_VM_Vector=zeros(M,n);
            
            r_opt_save=zeros(1,M);
            f_opt_save=zeros(1,M);
            L_opt_save=zeros(1,M);
            r_opt=r_opt_save;
            iter_consolidate_slot=zeros(1,n);
            NO_consolidate=zeros(1,n);
            
            
            %%%%%%%%%%%% values before consolidations
            E_W_old=zeros(nclient,n);                                                       %TCP/IP energy consumption for each slot
            Mu_old=zeros(nclient,n);  %[ones(nclient,1) zeros(nclient,n)];                                      %moltiplicatore
            E_TOT_old=zeros(1,n);                                                         %TOT energy consumption for each slot
            
            q_old=zeros(nclient,n);                                                        %service rate OUT-queue sistema
            s_old=zeros(nclient,n);                                                        %service rate IN-queue sistema
            a_old=zeros(nclient,n);                                                        %incoming rate IN-queue sistema
            r_old=zeros(nclient,n);                                                        %arrival rate in INPUT queue
            L_R_old=zeros(nclient,n);                                                      %real workload incoming to the buffer/queue
            L_tot_old=zeros(1,n);                                                          %workload per slot of the sistema
            lambda_old=zeros(nclient,n);                                                   %arrival rate in INPUT queue
            %L_tot_vect = 6+ round((10-6) * rand(1,n));
            
            %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
            %avg_r(nclient,1)=0.001;
            mu_save=0;
            for i=1:n                                                                   % SLOT
                
                if (mod(i,99999)==0)
                    i
                end
                
                %%%%%%% calculate avg_r
                if  (i>1)
                    
                    %%%%%%%%%%%%%%%%%NORMAL, without Jacobson
                    avg_r(nclient,i)= mean(r(nclient,1:1:i-1));                         %\overline(r)(t-1)
                    avg_r_old(nclient,i)= mean(r_old(nclient,1:1:i-1));                         %\overline(r)(t-1)
                    
                    %%%%%%%Jacobson formula
                    %                     prev_r=r(nclient,i-1);
                    %                     prev_r_old=r_old(nclient,i-1);
                    %                     avg_r(nclient,i)= (1-J_P)*avg_r(nclient,i-1)+J_P*prev_r;        %Jacoupson formula
                    %                     avg_r_old(nclient,i)= (1-J_P)*avg_r_old(nclient,i-1)+J_P*prev_r_old;        %Jacoupson formula
                end
                
                % calculate r (t)
                Mu_temp=(beta*(1-teta))./Mu(nclient,i);
                sigma_temp=((sigma(nclient,i).^2)./2).*Mu_temp;
                
                Mu_temp_old=(beta*(1-teta))./Mu_old(nclient,i);
                sigma_temp_old=((sigma(nclient,i).^2)./2).*Mu_temp_old;
                
                r(nclient,i)=min(max(sigma_temp,r_min),min(r_max,q(nclient,i)));
                r_old(nclient,i)=min(max(sigma_temp_old,r_min),min(r_max,q_old(nclient,i)));
                
                % calculate L_R (t) and buffer status
                %                 L_R_threshold=10^(-6);
                L_R(nclient,i)=max(r_min,avg_r(nclient,i));
                L_R_old(nclient,i)=max(r_min,avg_r_old(nclient,i));
                %if abs(r(nclient,i)-q(nclient,i))>L_R_threshold
                
                %OLD results (NO Consolidate)
                if r_old(nclient,i)==q_old(nclient,i)
                    empty_buffer_index1=1;
                else
                    empty_buffer_index1=0;
                end
                %NEW results (Consolidate)
                if r(nclient,i)==q(nclient,i)
                    empty_buffer_index=1;
                else
                    empty_buffer_index=0;
                end
                
                
                
                %%%%%%% calculate L_opt
                %OLD results (NO Consolidate)
                if (q_old(nclient,i)<=N_O+r_min-r_max)
                    
                    L_tot_old(i)= min(s_old(nclient,i),L_R_old(nclient,i)*(1+0.2*empty_buffer_index1));    %L_opt^(*)(t) with empty buffer state
                    
                    %                     L_tot_old(i)= min(s_old(nclient,i),max(r_min,avg_r_old(nclient,i)));             %L_opt^(*)(t)
                else
                    L_tot_old(i)= 0;
                end
                %NEW results (Consolidate)
                if (q(nclient,i)<=N_O+r_min-r_max)
                    L_tot(i)= min(s(nclient,i),L_R(nclient,i)*(1+0.2*empty_buffer_index));          %L_opt^(*)(t) with empty buffer state
                    
                    %                     L_tot(i)= min(s(nclient,i),max(r_min,avg_r(nclient,i)));                         %L_opt^(*)(t)
                    
                else
                    L_tot(i)= 0;
                    
                end
                
                % calculate lambda
                %XXXXXXXX
                lambda(nclient,i)=min(a(nclient,i),N_I-s(nclient,i)+L_tot(i));
                
                lambda_old(nclient,i)=min(a(nclient,i),N_I-s_old(nclient,i)+L_tot_old(i));
                
                
                
                % calculate E_W (t)
                if i>1
                    E_W(nclient,i)=(r(nclient,i).^2)./(sigma(nclient,i).^2);
                    E_W_old(nclient,i)=(r_old(nclient,i).^2)./(sigma(nclient,i).^2);
                else
                    E_W(nclient,1)=0;
                    E_W_old(nclient,1)=0;
                end
                % calculate RTT (t+1)
                % RTT(nclient,i) = (0.75 * RTT(nclient,i-1) + 0.25 * DELTA_IP1(nclient,i));
                
                % calculate Sigma (t+1)
                %  sigma(nclient,i) = ( K_0 * ((z(nclient,i))^0.5) ) / RTT(nclient,i);
                
                % calculate Mu(t+1)
                Mu(nclient,i+1)=max(0,Mu(nclient,i)+(kappa./i)*(E_W(nclient,i)-energy_ave));      %  min(limSupMu,max(0,Mu(nclient,i) + (energy_irradiata(nclient,i)-energy_ave)));
                Mu_old(nclient,i+1)=max(0,Mu_old(nclient,i)+(kappa./i)*(E_W_old(nclient,i)-energy_ave));      %  min(limSupMu,max(0,Mu(nclient,i) + (energy_irradiata(nclient,i)-energy_ave)));
                %Mu(nclient,i+1)=max(0,Mu(nclient,i)+(kappa./i)*(E_W(nclient,i)-energy_ave_vector(test1)));      %  min(limSupMu,max(0,Mu(nclient,i) + (energy_irradiata(nclient,i)-energy_ave)));
                % calculate q(t+1)
                q(nclient,i+1)=max(0,q(nclient,i)-r(nclient,i))+L_tot(i);
                q_old(nclient,i+1)=max(0,q_old(nclient,i)-r_old(nclient,i))+L_tot_old(i);
                
                % calculate s(t+1)
                s(nclient,i+1)=max(0,s(nclient,i)-L_tot(i))+lambda(nclient,i);
                s_old(nclient,i+1)=max(0,s_old(nclient,i)-L_tot_old(i))+lambda_old(nclient,i);
                
                % COMNET or DVFS approach
                % end % n
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %      for num_job=1:length(Jobs)
                % L_tot(i)=Jobs(num_job);
                
                %L_tot(i)=a(i);
                if (L_tot(i)==0)
                    mu_opt=0;
                    L_opt=zeros(1,M);
                    f_opt=f_zero;
                    iter=0;
                    
                    if (i==1)
                        active_VM(i)=0;
                        active_VM_old(i)=0;
                        active_VM_Vector(:,1)=Active_VM_rate;
                    else
                        active_VM(i)=active_VM(i-1);
                        active_VM_old(i)=active_VM_old(i-1);
                        active_VM_Vector(:,i)=active_VM_Vector(:,i-1);
                    end
                else
                    %CHECK FEASIBILITY - controllare funzionamento operatori %%%%%%%%%%%%%%%%%%
                    M1_vector(i)=M1;
                    if (i==1)
                        active_VM_Vector(:,1)=Active_VM_rate;
                        
                    else
                        active_VM_Vector(:,i)=active_VM_Vector(:,i-1);
                    end
                    condizione_feas_carico=(sum(f_max.*Delta-L_b)>=L_tot(i));
                    %condizione_feas_carico=(sum(active_VM_Vector(:,i)'.*f_max.*Delta-L_b)>=L_tot(i));
                    
                    condizione_feas_tempo=(L_tot(i)<=R_tot.*(T_tot-Delta)./2);
                    condizione_feas_back=min(Delta.*f_max>=L_b);
                    
                    
                    feasibility=condizione_feas_carico && condizione_feas_tempo && condizione_feas_back;
                    
                    if ~feasibility
                        % 'VM='
                        % M
                        SLA_violate(i)=1;  % the system unable to response the workload
                        active_VM(i)=active_VM(i-1);
                        active_VM_old(i)=active_VM_old(i-1);
                        active_VM_Vector(:,i)=active_VM_Vector(:,i-1);
                        W_range(i)=W_range(i-1);
                        W_range_old(i)=W_range(i);
                        
                        if (ceil(L_tot(i)/L_max(1)))>M % the workload is overload capacity and SLA violate and also estimate should considred
                            needed_new_VM(i)=ceil(L_tot(i)/(L_max(1)*up_th))-active_VM(i); % this is done just for hemonegnous case
                            L_opt_estimate=(L_tot(i)/(active_VM(i)+needed_new_VM(i))).*ones(1,active_VM(i)+needed_new_VM(i));
                            L_max_estimate=L_max(1).*ones(1,active_VM(i)+needed_new_VM(i));
                            W_range_estimate(i)=(1/(active_VM(i)+needed_new_VM(i)))*sum(L_opt_estimate./L_max_estimate);
                        end
                        continue
                        error('Problem unfeasible')
                    else
                        % 'Problem Feasibile!'
                    end
                    
                    
                    %   [conv,iter,mu_opt,ni_opt,f_opt,L_opt]=tracking(L_tot(i),L_b,alpha_zero,alpha_mu,chi,W,T_tot,Delta,f_zero,f_max,M,Th);
                    %[conv,iter,mu_opt,ni_opt,f_opt,L_opt,alpha_step,V_step]=tracking3(L_tot(i),L_b,alpha_zero,alpha_mu,chi,W,T_tot,Delta,f_zero,f_max,M,Th,mu_opt,ni_opt,f_opt,L_opt,alpha_step,V_step);
                    %[conv,iter,mu_opt,L_opt,alpha_step_finale,V_step_finale]= tracking_mu(L_tot(i),L_b,chi,W,T_tot,Delta,f_max,M,Th);
                    if (i==1)
                        [iter,mu_opt,delta_mu,L_opt]= Mu_opt_bisezione1(Active_VM_rate.*f_max,Delta,Active_VM_rate.*Th,L_tot(i),W,T_tot);
                    else
                        if (NO_consolidate(i-1)==0) % previous slot did not go to consolidation
                            if sum(Delta.*(Active_VM_rate.*f_max))>=L_tot(i) % IS system able to response the incoming workload?
                                [iter,mu_opt,delta_mu,L_opt]= Mu_opt_bisezione1(Active_VM_rate.*f_max,Delta,Active_VM_rate.*Th,L_tot(i),W,T_tot);
                            elseif (active_VM(i)>M)
                                SLA_violate(i)=1;  % the system unable to response the workload
                            elseif ((active_VM(i)<=M) && (ceil(L_tot(i)/L_max(1))>=M1)) % incoming workload is high but can respond by at most M VMS
                                new_added=ceil(L_tot(i)/L_max(1));
                                new_added=new_added-active_VM(i-1);
                                cc=0;
                                for ff=1:M
                                    if (Active_VM_rate(ff)==0) && cc<new_added
                                        Active_VM_rate(ff)=1;
                                        cc=cc+1;
                                    end
                                end
                                active_VM(i)=active_VM(i)+new_added;
                                [iter,mu_opt,delta_mu,L_opt]= Mu_opt_bisezione1(Active_VM_rate.*f_max,Delta,Active_VM_rate.*Th,L_tot(i),W,T_tot);
                            end
                        else
                            if (sum(Delta.*(Active_VM_rate.*f_max))>=L_tot(i)) && sum(Active_VM_rate)<=M
                                [iter,mu_opt,delta_mu,L_opt]= Mu_opt_bisezione1(Active_VM_rate.*f_max,Delta,Active_VM_rate.*Th,L_tot(i),W,T_tot);
                            elseif sum(Active_VM_rate)<=M
                                new_VM_updated=ceil(L_tot(i)/L_max(1));
                                new_VM_updated=new_VM_updated-sum(Active_VM_rate);
                                [row, index]=find(Active_VM_rate==0,new_VM_updated);
                                Active_VM_rate(index)=1;
                                [iter,mu_opt,delta_mu,L_opt]= Mu_opt_bisezione1(Active_VM_rate.*f_max,Delta,Active_VM_rate.*Th,L_tot(i),W,T_tot);
                            end
                        end
                    end
                    allocazione_parziale=L_opt<f_ibernazione.*Delta;
                    
                    %f_opt(L_opt<f_ibernazione.*Delta)=f_ibernazione(L_opt<f_ibernazione.*Delta);
                    %f_opt(L_opt>=f_ibernazione.*Delta)=L_opt(L_opt>=f_ibernazione.*Delta)./Delta;
                    f_opt(allocazione_parziale)=f_ibernazione(allocazione_parziale);
                    f_opt(~allocazione_parziale)=L_opt(~allocazione_parziale)./Delta;
                    
                    
                    
                    
                    
                    active_VM(i)=0;
                    active_VM(i)=sum(ne(f_opt,f_ibernazione));
                    active_VM_old(i)=active_VM(i);
                    active_VM_Vector(:,i)=ne(f_opt,f_ibernazione);
                    Active_VM_rate=(active_VM_Vector(:,i))';
                end % if
                
                W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                W_range_old(i)=W_range(i);
                
                down_th=1-up_th;
                if i>1
                    
                    VM_inactive_vector=~active_VM_Vector(:,i-1);
                    %VM_inactive_vector(:)=VM_inactive_vector(1:M1);
                    VM_inactive_vector(M1+1:M)=1;
                else
                    VM_inactive_vector=~Active_VM_rate';
                end
                if (W_range(i)> down_th) && (W_range(i)< up_th)
                    NO_consolidate(i)=0;  % no consolidation needed for non of servers/VMs
                    W_range_estimate(i)=W_range(i);
                    %elseif (W_range(i)> down_th) && (W_range(i)> up_th)
                    %   NO_consolidate(i)=2;  % ON servers/VMS
                    % [L_opt,f_opt,active_VM(i),L_b]=consolidate_type1(W_range(i),up_th, down_th,L_opt,L_opt_max,f_opt,f_ibernazione,Delta, active_VM(i),L_b);
                    %                         [consolidate_pos1, consolidate_value1]=max(L_opt);
                    %                         L_b(consolidate_value1)=L_opt(consolidate_value1);
                    %elseif (W_range(i)< down_th)
                elseif ((W_range(i)< down_th) || ((W_range(i)> up_th))) && (L_tot(i)>0)
                    NO_consolidate(i)=1; % OFF servers/VMS
                    
                    % [L_opt,f_opt,active_VM(i),L_b]=consolidate_type1(W_range(i),up_th, down_th,L_opt,L_opt_max,f_opt,f_ibernazione,Delta, active_VM(i),L_b);
                    %%consolidation phase
                    %                 L_opt_save=L_opt;
                    %                 mu_save=mu_opt;
                    %                 f_opt_save=f_opt;
                    %                 r_opt_save=r(nclient,i);
                    %
                    
                    gamma_max=0.4;
                    
                    TH=2*log(2)*N_o_i;
                    % NICOLA approach
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     [ iter_mu, iter_L_opt, iter_f_opt, L_opt, f_opt] = Consolidate_operation(f_opt_ON, r(nclient,i), f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                    %                      Delta, f_zero, VM_status_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save,mu_save, f_opt_save, r_opt_save, gamma_max);
                    %
                    %                     [iter,L_opt,delta_li,f_out]= L_opt_bisezione(L_max, M, N_o_i, W, T_tot, Delta, g, teta, alpha);
                    %                     [iter1,f_opt,f_f_opt,func_out]= f_opt_bisezione(f_opt_ON, r(nclient,i), f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                    %                     Delta, f_zero, VM_status_vector, M);
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % bacarelli approach
                    alpha=2;
                    if (i==1)
                        mu_opt1=0;
                        r_opt=zeros(1,M);
                    end
                    %                     [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli(f_opt_ON, f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                    %                         Delta, f_zero, VM_inactive_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH);
                    %
                    kappa1=100;
                    c1=10^(-4);
                    l_iteratiuons1=10^5;
                    tol_workload_allocated1=10^(-3);
                    %                     [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli1(f_opt_ON, f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                    %                         Delta, f_zero, VM_inactive_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH);
                    %
                    %                     [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(f_opt_ON, f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                    %                         Delta, f_zero, VM_inactive_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH,...
                    %                         kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                    %Active_VM_rate=ones(1,M);
                    if i>1 && (sum(NO_consolidate(1:i))==1)
                        Active_VM_rate=zeros(1,M);
                        %test=active_VM_Vector(:,i-1)';
                        Active_VM_rate(1:M1)=ones(1,M1);
                    elseif i>1
                        Active_VM_rate=zeros(1,M);
                        test=active_VM_Vector(:,i-1)';
                        Active_VM_rate(1:M1)=test(1:M1);
                    end
                    [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(Active_VM_rate.*f_opt_ON, Active_VM_rate.*f_max, k_e, teta, alpha, g, sigma1, Active_VM_rate.*E_idle, Active_VM_rate.*E_max, T_ON,...
                        Delta, Active_VM_rate.*f_zero, VM_inactive_vector, M, Active_VM_rate.*L_max, L_tot(i), N_o_i, T_tot,W, Active_VM_rate'.*L_opt_save',mu_save, Active_VM_rate'.*f_opt_save', r_opt_save, gamma_max, TH,...
                        kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                    
                    
                    %[iter,li,delta_li,f_out]= L_opt_bisezione(L_max,M, N_o_i, W, T_tot, Delta, g, teta, alpha, mu_opt1);
                    iter_consolidate_slot(i)=iter_consolidate;
                    f_opt=f_opt';
                    L_opt=L_opt';
                    %r(nclient,i)=min(r_min,r_opt(1);
                    m_opt=mu_opt1;
                    %L_opt=li;
                    
                    if ceil(L_tot(i)/L_max(1))<M  && sum(L_max(1:ceil(L_tot(i)/L_max(1))))>=L_tot(i) % just for Hemogenous
                        L_opt=zeros(1,M);
                        L_opt(1:ceil(L_tot(i)/L_max(1)))=L_max(1:ceil(L_tot(i)/L_max(1)));
                    end
                    
                    f_opt=L_opt./Delta;
                    
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    L_ibernazione=f_ibernazione*Delta;
                    active_VM_Vector(:,i)=ne(f_opt,f_ibernazione)& ne(L_opt,L_ibernazione);
                    active_VM(i)=sum(active_VM_Vector(:,i));
                    W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                    Active_VM_rate=active_VM_Vector(:,i)';
                    % M1=active_VM(i);
                    %                     temp_active=sum(~VM_inactive_vector);
                    %                     temp_active_vector=~VM_inactive_vector;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % this is done just for hemonegnous case
                    
                    % we are in high range and uanble to decrease W_rang
                    if ((W_range(i)>= up_th) && (W_range(i)>W_range_old(i)) && active_VM(i)<M && M1<M)
                        while (W_range(i)>= up_th) && (W_range(i)>W_range_old(i)) && W_range(i)<1 && active_VM(i)<=M
                            if (active_VM(i)<=M1) && ((active_VM(i)<M))
                                active_VM(i)=active_VM(i)+1;
                                L_opt=zeros(1,M);
                                L_opt(1:active_VM(i))=(L_tot(i)/active_VM(i)).*ones(1,active_VM(i));
                                f_opt=L_opt./Delta;
                                W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                                M1=active_VM(i);
                                active_VM_Vector(1:M1,i)=ones(1,M1);
                                Active_VM_rate=active_VM_Vector(:,i)';
                                
                            end
                        end
                    end
                    % we are in high range and uanble to decrease W_rang
                    if ((W_range(i)>= up_th) && (W_range(i)<W_range_old(i)) && active_VM(i)<M && M1<M)
                        while (W_range(i)>= up_th) && (W_range(i)<W_range_old(i) && (W_range(i)>down_th) && active_VM(i)<=M)
                            if (active_VM(i)<=M1) && ((active_VM(i)<M))
                                active_VM(i)=active_VM(i)+1;
                                L_opt=zeros(1,M);
                                L_opt(1:active_VM(i))=(L_tot(i)/active_VM(i)).*ones(1,active_VM(i));
                                f_opt=L_opt./Delta;
                                W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                                M1=active_VM(i);
                                active_VM_Vector(1:M1,i)=ones(1,M1);
                                Active_VM_rate=active_VM_Vector(:,i)';
                            end
                            if (active_VM(i)==M) % we are in high range and unable to go out after decreament and VM is full working
                                break
                            end
                        end
                    end
                    %estimation
                    if ((W_range(i)>= up_th) && (W_range(i)<W_range_old(i)))...
                            || ((W_range(i)>= up_th) && (W_range(i)==W_range_old(i)))% we are in high range and decrease W_rang but unable to go to the normal environment
                        % this is done just for hemonegnous case
                        needed_new_VM(i)=ceil(L_tot(i)/(L_max(1)*up_th))-active_VM(i); % this is done just for hemonegnous case
                        
                        %L_opt_estimate=(L_tot(i)/(active_VM(i)+needed_new_VM(i))).*ones(1,M+needed_new_VM(i));%old
                        %style
                        %
                        
                        L_opt_estimate=(L_tot(i)/(active_VM(i)+needed_new_VM(i))).*ones(1,M+needed_new_VM(i));
                        L_max_estimate=L_max(1).*ones(1,M+needed_new_VM(i));
                        W_range_estimate(i)=(1/(active_VM(i)+needed_new_VM(i)))*sum(L_opt_estimate./L_max_estimate);
                    elseif ((W_range(i)>= up_th) && (W_range(i)==1) && active_VM(i)<M)
                        while ((W_range(i)>= up_th)&& active_VM(i)<M) % we have enough OFF VMs and we can turn on to decrease W_range
                            active_VM(i)=active_VM(i)+1;
                            L_opt=zeros(1,M);
                            %                             kappa1=1000;
                            %                             alpha=2;
                            %                             c1=10^(-4);
                            %                             l_iteratiuons1=10^5;
                            %                             tol_workload_allocated1=10^(-3);
                            %                             [row, index]=find(VM_inactive_vector==1,1);
                            %
                            %                             VM_inactive_vector(row)=0;
                            %                                 [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(f_opt_ON, (~VM_inactive_vector)'.*f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                            %                                     Delta, f_zero, VM_inactive_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH,...
                            %                                     kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                            %                                 if (iter_consolidate<l_iteratiuons1)
                            %                                 [row, index]=find(VM_inactive_vector==1,3);
                            %                                 VM_inactive_vector(row)=0;
                            %                                 [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(f_opt_ON, (~VM_inactive_vector)'.*f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                            %                                     Delta, f_zero, VM_inactive_vector, M, L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH,...
                            %                                     kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                            %                                 end
                            %                                 iter_consolidate_slot(i)=iter_consolidate;
                            %                                 f_opt=f_opt';
                            %                                 L_opt=L_opt';
                            %                                 r(nclient,i)=r_opt(1);
                            %                                 m_opt=mu_opt1;
                            %                                 %L_opt=li;
                            %                                 f_opt=L_opt./Delta;
                            %                                 L_ibernazione=f_ibernazione*Delta;
                            %                                 active_VM_Vector(:,i)=ne(f_opt,f_ibernazione)& ne(L_opt,L_ibernazione);
                            %                                 active_VM(i)=sum(active_VM_Vector(:,i));
                            
                            
                            L_opt(1:active_VM(i))=(L_tot(i)/active_VM(i)).*ones(1,active_VM(i));
                            f_opt=L_opt./Delta;
                            W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                            M1=active_VM(i);
                            active_VM_Vector(1:M1,i)=ones(1,M1);
                            Active_VM_rate=active_VM_Vector(:,i)';
                            W_range_estimate(i)=W_range(i); % for calculating the estimation to reach the up-range
                        end
                        % estimation to lower the W. VM is the maximum and
                        % L_tot is too much and unable to reach the normal
                        % range
                        if ((W_range(i)>= up_th)&& active_VM(i)==M)
                            needed_new_VM(i)=ceil(L_tot(i)/(L_max(1)*up_th))-active_VM(i); % this is done just for hemonegnous case
                            L_opt_estimate=(L_tot(i)/(active_VM(i)+needed_new_VM(i))).*ones(1,M+needed_new_VM(i));
                            L_max_estimate=L_max(1).*ones(1,M+needed_new_VM(i));
                            W_range_estimate(i)=(1/(active_VM(i)+needed_new_VM(i)))*sum(L_opt_estimate./L_max_estimate);
                        end
                    else
                        W_range_estimate(i)=W_range(i); % for calculating the estimation to reach the up-range
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%checking
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for turn
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OFF
                    tt=1;zz=2;
                    while_step=0;
                    while ((W_range(i)<= down_th) && (W_range(i)<W_range_old(i))) || ((W_range(i)<= down_th) && (while_step<M)) ...
                            || (up_th <= W_range(i)&& (W_range(i) < 1) && (W_range(i)>W_range_old(i)))
                        while_step=while_step+1;
                        if (W_range(i)>down_th) &&(W_range(i) < 1) && (W_range(i)>W_range_old(i)) && (while_step<M) % go out from the lower case and come to the normal case
                            break
                        end
                        if ((W_range(i)<= down_th) && (W_range(i)<W_range_old(i))) % still stay lower that low range
                            kappa1=1000;
                            alpha=2;
                            c1=10^(-4);
                            l_iteratiuons1=10^5;
                            tol_workload_allocated1=10^(-3);
                            [row, index]=find(VM_inactive_vector==0,tt);
                            
                            VM_inactive_vector(row)=1;
                        elseif ((W_range(i)<= down_th) && (W_range(i)>W_range_old(i))) % still stay lower that low range
                            kappa1=10000;
                            alpha=2;
                            c1=10^(-4);
                            l_iteratiuons1=10^5;
                            tol_workload_allocated1=10^(-3);
                            [row, index]=find(VM_inactive_vector==0,tt);
                            VM_inactive_vector(row)=1;
                        elseif ((W_range(i)> up_th)&& (W_range(i)< 1)) % still stay higher than high range
                            kappa1=1000;
                            alpha=1;
                            c1=10^(-2);
                            l_iteratiuons1=10^5;
                            tol_workload_allocated1=10^(-3);
                            %                             for hi=1:M
                            %                                 if (sum(f_max(1:hi).*Delta)>=L_tot(i)) % find the most
                            %                                     break
                            %                                 end
                            %                             end
                            
                            
                        end
                        
                        
                        %active_VM_Vector(M-tt:M,i)= 0;
                        %VM_inactive_vector=~active_VM_Vector(:,i);
                        
                        if sum(VM_inactive_vector)==M  % if consolidation turned whol system off we should avoid and put the first ON
                            VM_inactive_vector(1)=0;
                            %zz=zz+1;
                        end
                        %                         if (W_range(i)> up_th)
                        %                             VM_inactive_vector(zz)=0;
                        %                         end
                        
                        [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(f_opt_ON, (~VM_inactive_vector)'.*f_max, k_e, teta, alpha, g, sigma1, E_idle, E_max, T_ON,...
                            Delta, (~VM_inactive_vector)'.*f_zero, VM_inactive_vector, M, (~VM_inactive_vector)'.*L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH,...
                            kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                        %                 [ iter_consolidate, L_opt, f_opt, r_opt, mu_opt1] = Consolidate_operation_baccarelli2(f_opt_ON, Active_VM_rate.*f_max, k_e, teta, alpha, g, sigma1, Active_VM_rate.*E_idle, Active_VM_rate.*E_max, T_ON,...
                        %                     Delta, Active_VM_rate.*f_zero, VM_inactive_vector, M, Active_VM_rate.*L_max, L_tot(i), N_o_i, T_tot,W, L_opt_save',mu_save, f_opt_save', r_opt_save, gamma_max, TH,...
                        %                     kappa1, c1, tol_workload_allocated1, l_iteratiuons1);
                        
                        iter_consolidate_slot(i)=iter_consolidate;
                        f_opt=f_opt';
                        L_opt=L_opt';
                        %r(nclient,i)=r_opt(1);
                        m_opt=mu_opt1;
                        %L_opt=li;
                        
                        
                        f_opt=L_opt./Delta;
                        L_ibernazione=f_ibernazione*Delta;
                        active_VM_Vector(:,i)=ne(f_opt,f_ibernazione)& ne(L_opt,L_ibernazione);
                        active_VM(i)=sum(active_VM_Vector(:,i));
                        W_range(i)=(1/active_VM(i))*sum(L_opt./L_max);
                        W_range_estimate(i)=W_range(i); % for calculating the estimation to reach the up-range
                        M1_new=active_VM(i);
                        M1=M1_new;
                        Active_VM_rate=active_VM_Vector(:,i)';
                        
                        tt=tt+1;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%checking
                    end
                end
                
                
                %%
                
                
                R_opt=2.*L_opt./(T_tot-Delta);
                iter_VM(k)=iter_VM(k)+iter;
                
                
                
                L_VM_prima(k)=L_VM_prima(k)+L_opt(1);
                L_VM_ultima(k)=L_VM_ultima(k)+L_opt(M);
                f_VM_prima(k)=f_VM_prima(k)+f_opt(1);
                f_VM_ultima(k)=f_VM_ultima(k)+f_opt(M);
                
                
                
                %Calcolo f_precedente e f_successiva a partire da
                %f_opt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                f_precedente=zeros(1,M);
                f_successiva=zeros(1,M);
                x=zeros(1,M);
                
                for conta=1:M
                    
                    delta_f_discrete=f_discrete(conta,:)-f_opt(conta);
                    [ff,ind_ff]=min(abs(delta_f_discrete));
                    if ff==0
                        f_precedente(conta)=f_discrete(conta,ind_ff);
                        f_successiva(conta)=f_discrete(conta,ind_ff);
                        x(conta)=1; %qualsiasi valore è indifferente
                    elseif ind_ff==1
                        f_precedente(conta)=f_discrete(conta,1);
                        f_successiva(conta)=f_discrete(conta,2);
                        x(conta)=(f_opt(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
                    elseif ind_ff==Q
                        f_precedente(conta)=f_discrete(conta,Q-1);
                        f_successiva(conta)=f_discrete(conta,Q);
                        x(conta)=(f_opt(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
                    elseif delta_f_discrete(ind_ff)>0
                        f_precedente(conta)=f_discrete(conta,ind_ff-1);
                        f_successiva(conta)=f_discrete(conta,ind_ff);
                        x(conta)=(f_opt(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
                    else
                        f_precedente(conta)=f_discrete(conta,ind_ff);
                        f_successiva(conta)=f_discrete(conta,ind_ff+1);
                        x(conta)=(f_opt(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%  costo  %%%%%%%%%%%%%%%%%%
                %                     if (~conv)
                %                         error('Errore: Carico non correttamente allocato')
                %                     end
                
                costo_speed=sum(((f_opt./f_max).^2).*omega.*E_max);
                costo_switch=sum(k_e.*(f_opt-f_zero).^2);
                %costo_channel=sum(2.*P_net.*L_opt./C_max);
                costo_channel=(T_tot-Delta).*sum(chi.*((2.^(R_opt./W))-1));
                %costo_tot=costo_speed+costo_switch+costo_channel+E_W(i);
                costo_tot=costo_speed+costo_switch+costo_channel;
                
                
                
                %%%%%%%%%%%%%%%%%%   COSTI FREQUENZE DISCRETE  %%%%%%%%%
                
                
                % Costo frequenza_successiva
                costo_speed_f_succ=sum(((f_successiva./f_max).^2).*omega.*E_max);
                costo_switch_f_succ=sum(k_e.*(f_successiva-f_zero).^2);
                costo_channel_f_succ=(T_tot-Delta).*sum(chi.*((2.^(R_opt./W))-1));
                %costo_tot_f_prec=costo_speed_f_succ+costo_switch_f_succ+costo_channel_f_succ+E_W(i);
                costo_tot_f_succ=costo_speed_f_succ+costo_switch_f_succ+costo_channel_f_succ;
                
                % Costo frequenza_precedente
                costo_speed_f_prec=sum(((f_precedente./f_max).^2).*omega.*E_max);
                costo_switch_f_prec=sum(k_e.*(f_precedente-f_zero).^2);
                costo_channel_f_prec=(T_tot-Delta).*sum(chi.*((2.^(R_opt./W))-1));
                %costo_tot_f_prec=costo_speed_f_prec+costo_switch_f_prec+costo_channel_f_prec+E_W(i);
                costo_tot_f_prec=costo_speed_f_prec+costo_switch_f_prec+costo_channel_f_prec;
                
                % Costo Time-sharing
                P_Idle=10.*ones(1,active_VM(i)); %watt
                
                costo_speed_TS=sum((((f_precedente./f_max).^2).*omega.*E_max).*(1-x)+(((f_successiva./f_max).^2).*omega.*E_max).*x)+sum(P_Idle);
                costo_switch_TS=sum((k_e.*(f_precedente-f_zero).^2).*(1-x)+(k_e.*(f_successiva-f_zero).^2).*x);
                costo_channel_TS=(T_tot-Delta).*sum(chi.*((2.^(R_opt./W))-1));
                %costo_tot_TS=costo_speed_TS+costo_switch_TS+costo_channel_TS+E_W(i);
                costo_tot_TS=costo_speed_TS+costo_switch_TS+costo_channel_TS;
                
                
                %Costo Time-sharing-switch for heterogenous case
                %costo_TS_switch=sum(k_e.*(f_successiva-f_precedente).^2);
                %costo_tot_TSpiuSwitch1=costo_tot_TS+costo_TS_switch;
                costo_tot_TSpiuSwitch1=costo_tot_TS;
                
                
                
                %%%%%%%%%%%%%%%%%%   FINE  COSTI FREQUENZE DISCRETE  %%%%%%%%%
                
                
                
                costo_tot_VM(k)=costo_tot_VM(k)+costo_tot;
                costo_tot_VM_f_succ(k)=costo_tot_VM_f_succ(k)+costo_tot_f_prec;
                costo_tot_VM_f_prec(k)=costo_tot_VM_f_prec(k)+costo_tot_f_prec;
                costo_tot_VM_TS(k)=costo_tot_VM_TS(k)+costo_tot_TS;
                %costo_VM_TS_switch(k)=costo_VM_TS_switch(k)+costo_TS_switch;
                costo_tot_VM_TSpiuSwitch(k)=costo_tot_VM_TSpiuSwitch(k)+costo_tot_TSpiuSwitch1;
                %E_TOT(nclient,i)=costo_tot_TSpiuSwitch1;
                E_TOT(i)=costo_tot_TSpiuSwitch1;
                
                
                
                L_opt_save=L_opt;
                mu_save=mu_opt;
                f_opt_save=f_opt;
                %r_opt_save=r(nclient,i);
                r_opt_save=r_opt(1);
                
                f_zero=f_opt;
                
                
            end       %for each slot
            
            %         E_W_vector(k,test1,test2,test3,test31)=mean(E_W(nclient,:));
            %         E_TOT_vector(k,test1,test2,test3,test31)=mean(E_TOT(nclient,:));
            %         E_TOT=zeros(nclient,n);
            %         Muinf_vector(k,test1,test2,test3,test31)=Mu(nclient,n+1);
            %         s_vector(k,test1,test2,test3,test31)=mean(s(nclient,:));
            %         q_vector(k,test1,test2,test3,test31)=mean(q(nclient,:));
            %         r_vector(k,test1,test2,test3,test31)=mean(r(nclient,:));
            %         L_tot_vector(k,test1,test2,test3,test31)=mean(L_tot(nclient,:));
            %         lambda_vector(k,test1,test2,test3,test31)=mean(lambda(nclient,:));
            %         Utility_vector(k,test1,test2,test3,test31)=((test31.*E_TOT_vector(k,test1,test2,test3,test31))-(beta.*(1-test31).*r_vector(k,test1,test2,test3,test31)));
            
            
            %             E_W_vector(test2, test11)=mean(E_W(nclient,:));
            %             E_TOT_vector(test2, test11)=mean(E_TOT(:));
            %             E_TOT=zeros(1,n);
            %             Muinf_vector(test2, test11)=Mu(nclient,n+1);
            %             s_vector(test2, test11)=mean(s(nclient,:));
            %             q_vector(test2, test11)=mean(q(nclient,:));
            %             r_vector(test2, test11)=mean(r(nclient,:));
            %             L_tot_vector(test2, test11)=mean(L_tot(nclient,:));
            %             lambda_vector(test2, test11)=mean(lambda(nclient,:));
            %             Utility_vector(test2, test11)=((teta.*E_TOT_vector(test2, test11))-(beta.*(1-teta).*r_vector(test2, test11)));
            %
            %         E_W_vector(test2)=mean(E_W(nclient,:));
            %         E_TOT_vector(test2)=mean(E_TOT(:));
            %         E_TOT=zeros(1,n);
            %         Muinf_vector(test2)=Mu(nclient,n+1);
            %         s_vector(test2)=mean(s(nclient,:));
            %         q_vector(test2)=mean(q(nclient,:));
            %         r_vector(test2)=mean(r(nclient,:));
            %         L_tot_vector(test2)=mean(L_tot(nclient,:));
            %         lambda_vector(test2)=mean(lambda(nclient,:));
            %           Utility_vector(test2, test11)=((teta.*E_TOT_vector(test2, test11))-(beta.*(1-teta).*r_vector(test2, test11)));
        end  % test11 J_P
    end % test2 various q0
    
    %end % test1 various sigma0
    
    
    %end % test1 various E_ave
    
    
    
    % E_TOT_vector(k)=mean(E_TOT(:));
    costo_tot_VM(k)=costo_tot_VM(k)./length(Jobs);
    costo_tot_VM_f_succ(k)=costo_tot_VM_f_succ(k)./length(Jobs);
    costo_tot_VM_f_prec(k)=costo_tot_VM_f_prec(k)./length(Jobs);
    costo_tot_VM_TS(k)=costo_tot_VM_TS(k)./length(Jobs);
    %costo_VM_TS_switch(k)=costo_VM_TS_switch(k)./length(Jobs);
    costo_tot_VM_TSpiuSwitch(k)=costo_tot_VM_TSpiuSwitch(k)./length(Jobs);
    
    
    iter_VM(k)=iter_VM(k)./length(Jobs);
    
    L_VM_prima(k)=L_VM_prima(k)./length(Jobs);
    L_VM_ultima(k)=L_VM_ultima(k)./length(Jobs);
    f_VM_prima(k)=f_VM_prima(k)./length(Jobs);
    f_VM_ultima(k)=f_VM_ultima(k)./length(Jobs);
    
    clear f_max f_zero Th E_max omega Delta L_b temp alpha_zero alpha_mu chi W ni_opt f_opt L_opt mu_opt alpha_step V_step f_precedente f_successiva delta_f_discrete x;
    
    toc
    
    
end    % various VM
%
%% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%testing
% r for 5 selected VMS types, for E_ave=0.75
% x11=[r_vector(1,1,1,1,1),r_vector(1,1,1,1,2),r_vector(1,1,1,1,3),r_vector(1,1,1,1,4),r_vector(1,1,1,1,5)];% VM=2
% x12=[r_vector(4,1,1,1,1),r_vector(4,1,1,1,2),r_vector(4,1,1,1,3),r_vector(4,1,1,1,4),r_vector(4,1,1,1,5)];% VM=5
% x13=[r_vector(9,1,1,1,1),r_vector(9,1,1,1,2),r_vector(9,1,1,1,3),r_vector(9,1,1,1,4),r_vector(9,1,1,1,5)];%VM=10
% x14=[r_vector(14,1,1,1,1),r_vector(14,1,1,1,2),r_vector(14,1,1,1,3),r_vector(14,1,1,1,4),r_vector(14,1,1,1,5)];%VM=15
% x15=[r_vector(19,1,1,1,1),r_vector(19,1,1,1,2),r_vector(19,1,1,1,3),r_vector(19,1,1,1,4),r_vector(19,1,1,1,5)];%VM=20
% r for 5 selected VMS types, for E_ave=0.5
% x21=[r_vector(1,2,1,1,1),r_vector(1,2,1,1,2),r_vector(1,2,1,1,3),r_vector(1,2,1,1,4),r_vector(1,2,1,1,5)];% VM=2
% x22=[r_vector(4,2,1,1,1),r_vector(4,2,1,1,2),r_vector(4,2,1,1,3),r_vector(4,2,1,1,4),r_vector(4,2,1,1,5)];% VM=5
% x23=[r_vector(9,2,1,1,1),r_vector(9,2,1,1,2),r_vector(9,2,1,1,3),r_vector(9,2,1,1,4),r_vector(9,2,1,1,5)];%VM=10
% x24=[r_vector(14,2,1,1,1),r_vector(14,2,1,1,2),r_vector(14,2,1,1,3),r_vector(14,2,1,1,4),r_vector(14,2,1,1,5)];%VM=15
% x25=[r_vector(19,2,1,1,1),r_vector(19,2,1,1,2),r_vector(19,2,1,1,3),r_vector(19,2,1,1,4),r_vector(19,2,1,1,5)];%VM=20
%
%
% r for 5 selected VMS types, for E_ave=0.25
% x31=[r_vector(1,3,1,1,1),r_vector(1,3,1,1,2),r_vector(1,3,1,1,3),r_vector(1,3,1,1,4),r_vector(1,3,1,1,5)];% VM=2
% x32=[r_vector(4,3,1,1,1),r_vector(4,3,1,1,2),r_vector(4,3,1,1,3),r_vector(4,3,1,1,4),r_vector(4,3,1,1,5)];% VM=5
% x33=[r_vector(9,3,1,1,1),r_vector(9,3,1,1,2),r_vector(9,3,1,1,3),r_vector(9,3,1,1,4),r_vector(9,3,1,1,5)];%VM=10
% x34=[r_vector(14,3,1,1,1),r_vector(14,3,1,1,2),r_vector(14,3,1,1,3),r_vector(14,3,1,1,4),r_vector(14,3,1,1,5)];%VM=15
% x35=[r_vector(19,3,1,1,1),r_vector(19,3,1,1,2),r_vector(19,3,1,1,3),r_vector(19,3,1,1,4),r_vector(19,3,1,1,5)];%VM=20
%
% r for 5 selected VMS types, for E_ave=0.125
% x41=[r_vector(1,4,1,1,1),r_vector(1,4,1,1,2),r_vector(1,4,1,1,3),r_vector(1,4,1,1,4),r_vector(1,4,1,1,5)];% VM=2
% x42=[r_vector(4,4,1,1,1),r_vector(4,4,1,1,2),r_vector(4,4,1,1,3),r_vector(4,4,1,1,4),r_vector(4,4,1,1,5)];% VM=5
% x43=[r_vector(9,4,1,1,1),r_vector(9,4,1,1,2),r_vector(9,4,1,1,3),r_vector(9,4,1,1,4),r_vector(9,4,1,1,5)];%VM=10
% x44=[r_vector(14,4,1,1,1),r_vector(14,4,1,1,2),r_vector(14,4,1,1,3),r_vector(14,4,1,1,4),r_vector(14,4,1,1,5)];%VM=15
% x45=[r_vector(19,4,1,1,1),r_vector(19,4,1,1,2),r_vector(19,4,1,1,3),r_vector(19,4,1,1,4),r_vector(19,4,1,1,5)];%VM=20
%
%
% E_tot for 5 selected VMS types, for E_ave=20
% y11=[E_TOT_vector(1,1,1,1,1),E_TOT_vector(1,1,1,1,2),E_TOT_vector(1,1,1,1,3),E_TOT_vector(1,1,1,1,4),E_TOT_vector(1,1,1,1,5)];% VM=2
% y12=[E_TOT_vector(4,1,1,1,1),E_TOT_vector(4,1,1,1,2),E_TOT_vector(4,1,1,1,3),E_TOT_vector(4,1,1,1,4),E_TOT_vector(4,1,1,1,5)];% VM=5
% y13=[E_TOT_vector(2,1,1,1,1),E_TOT_vector(2,1,1,1,2),E_TOT_vector(2,1,1,1,3),E_TOT_vector(2,1,1,1,4),E_TOT_vector(2,1,1,1,5)];%VM=10
% y14=[E_TOT_vector(14,1,1,1,1),E_TOT_vector(14,1,1,1,2),E_TOT_vector(14,1,1,1,3),E_TOT_vector(14,1,1,1,4),E_TOT_vector(14,1,1,1,5)];%VM=15
% y15=[E_TOT_vector(12,1,1,1,1),E_TOT_vector(12,1,1,1,2),E_TOT_vector(12,1,1,1,3),E_TOT_vector(12,1,1,1,4),E_TOT_vector(12,1,1,1,5)];%VM=20
% E_tot for 5 selected VMS types, for E_ave=10
% y21=[E_TOT_vector(1,2,1,1,1),E_TOT_vector(1,2,1,1,2),E_TOT_vector(1,2,1,1,3),E_TOT_vector(1,2,1,1,4),E_TOT_vector(1,2,1,1,5)];% VM=2
% y22=[E_TOT_vector(4,2,1,1,1),E_TOT_vector(4,2,1,1,2),E_TOT_vector(4,2,1,1,3),E_TOT_vector(4,2,1,1,4),E_TOT_vector(4,2,1,1,5)];% VM=5
% y23=[E_TOT_vector(9,2,1,1,1),E_TOT_vector(9,2,1,1,2),E_TOT_vector(9,2,1,1,3),E_TOT_vector(9,2,1,1,4),E_TOT_vector(9,2,1,1,5)];%VM=10
% y24=[E_TOT_vector(14,2,1,1,1),E_TOT_vector(14,2,1,1,2),E_TOT_vector(14,2,1,1,3),E_TOT_vector(14,2,1,1,4),E_TOT_vector(14,2,1,1,5)];%VM=15
% y25=[E_TOT_vector(19,2,1,1,1),E_TOT_vector(19,2,1,1,2),E_TOT_vector(19,2,1,1,3),E_TOT_vector(19,2,1,1,4),E_TOT_vector(19,2,1,1,5)];%VM=20
%
%
% E_tot for 5 selected VMS types, for E_ave=1
% y31=[E_TOT_vector(1,3,1,1,1),E_TOT_vector(1,3,1,1,2),E_TOT_vector(1,3,1,1,3),E_TOT_vector(1,3,1,1,4),E_TOT_vector(1,3,1,1,5)];% VM=2
% y32=[E_TOT_vector(4,3,1,1,1),E_TOT_vector(4,3,1,1,2),E_TOT_vector(4,3,1,1,3),E_TOT_vector(4,3,1,1,4),E_TOT_vector(4,3,1,1,5)];% VM=5
% y33=[E_TOT_vector(9,3,1,1,1),E_TOT_vector(9,3,1,1,2),E_TOT_vector(9,3,1,1,3),E_TOT_vector(9,3,1,1,4),E_TOT_vector(9,3,1,1,5)];%VM=10
% y34=[E_TOT_vector(14,3,1,1,1),E_TOT_vector(14,3,1,1,2),E_TOT_vector(14,3,1,1,3),E_TOT_vector(14,3,1,1,4),E_TOT_vector(14,3,1,1,5)];%VM=15
% y35=[E_TOT_vector(19,3,1,1,1),E_TOT_vector(19,3,1,1,2),E_TOT_vector(19,3,1,1,3),E_TOT_vector(19,3,1,1,4),E_TOT_vector(19,3,1,1,5)];%VM=20
%
% E_tot for 5 selected VMS types, for E_ave=0.75
% y41=[E_TOT_vector(1,4,1,1,1),E_TOT_vector(1,4,1,1,2),E_TOT_vector(1,4,1,1,3),E_TOT_vector(1,4,1,1,4),E_TOT_vector(1,4,1,1,5)];% VM=2
% y42=[E_TOT_vector(4,4,1,1,1),E_TOT_vector(4,4,1,1,2),E_TOT_vector(4,4,1,1,3),E_TOT_vector(4,4,1,1,4),E_TOT_vector(4,4,1,1,5)];% VM=5
% y43=[E_TOT_vector(9,4,1,1,1),E_TOT_vector(9,4,1,1,2),E_TOT_vector(9,4,1,1,3),E_TOT_vector(9,4,1,1,4),E_TOT_vector(9,4,1,1,5)];%VM=10
% y44=[E_TOT_vector(14,4,1,1,1),E_TOT_vector(14,4,1,1,2),E_TOT_vector(14,4,1,1,3),E_TOT_vector(14,4,1,1,4),E_TOT_vector(14,4,1,1,5)];%VM=15
% y45=[E_TOT_vector(19,4,1,1,1),E_TOT_vector(19,4,1,1,2),E_TOT_vector(19,4,1,1,3),E_TOT_vector(19,4,1,1,4),E_TOT_vector(19,4,1,1,5)];%VM=20
%
%
% figure(12000)
% subplot(2,2,1);
% plot(x11,y11,'--+',x12,y12,'--*',x13,y13,'--o',x14,y14,'--+',x15,y15,'--*');
% xlabel('$\overline{r}$ $(byte/slot)$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E}_{TOT}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('M=2','M=5','M=10','M=15','M=20');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $E_{ave}=0.75$, $\theta=\{0.01, 0.1, 0.25,0.5, 0.75, 0.9, 0.99\}$','Interpreter','latex','FontSize',12);
% grid on
%
% subplot(2,2,2);
% plot(x21,y21,'--+',x22,y22,'--*',x23,y23,'--o',x24,y24,'--+',x25,y25,'--*');
% xlabel('$\overline{r}$ $(byte/slot)$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E}_{TOT}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('M=2','M=5','M=10','M=15','M=20');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $E_{ave}=0.5$, $\theta=\{0.01, 0.1, 0.25,0.5, 0.75, 0.9, 0.99\}$','Interpreter','latex','FontSize',12);
% grid on
%
% subplot(2,2,3);
% plot(x31,y31,'--+',x32,y32,'--*',x33,y33,'--o',x34,y34,'--+',x35,y35,'--*');
% xlabel('$\overline{r}$ $(byte/slot)$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E}_{TOT}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('M=2','M=5','M=10','M=15','M=20');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $E_{ave}=0.25$, $\theta=\{0.01, 0.1, 0.25,0.5, 0.75, 0.9, 0.99\}$','Interpreter','latex','FontSize',12);
% grid on
%
% subplot(2,2,4);
% plot(x41,y41,'--+',x42,y42,'--*',x43,y43,'--o',x44,y44,'--+',x45,y45,'--*');
% xlabel('$\overline{r}$ $(byte/slot)$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E}_{TOT}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('M=2','M=5','M=10','M=15','M=20');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $E_{ave}=0.125$, $\theta=\{0.01, 0.1, 0.25,0.5, 0.75, 0.9, 0.99\}$','Interpreter','latex','FontSize',12);
% grid on
%
% figure(1)
% plot(VM,E_TOT_vector(:),'r');
% xlabel('$VMs$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_{TOT}}$ (Joule)','Interpreter','latex','FontSize',14);
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7000)
% subplot(2,2,1);
% plot(r_min_vector,E_TOT_vector(:),'r');
% xlabel('$r_{min}$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_{TOT}}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.125, \theta=0.5');
% title('$q_0=N_O=240$, $s_{0}=240$, $\sigma_{min}=4.8$, $r^{anal}_{min}=1.71$, $r_{min}=\sigma_{min}\sqrt{E_{ave}}$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(r_min_vector,E_W_vector(:),'r');
% xlabel('$r_{min}$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_{W}}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.125, \theta=0.5');
% title('$q_0=N_O=240$, $s_{0}=240$, $\sigma_{min}=4.8$, $r^{anal}_{min}=1.71$, $r_{min}=\sigma_{min}\sqrt{E_{ave}}$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(r_min_vector,r_vector(:),'r');
% xlabel('$r_{min}$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.125, \theta=0.5');
% title('$q_0=N_O=240$, $s_{0}=240$, $\sigma_{min}=4.8$, $r^{anal}_{min}=1.71$, $r_{min}=\sigma_{min}\sqrt{E_{ave}}$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(r_min_vector,Utility_vector(:),'r');
% xlabel('$r_{min}$','Interpreter','latex','FontSize',14);
% ylabel('$Utility$','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.125, \theta=0.5');
% title('$q_0=N_O=240$, $s_{0}=240$, $\sigma_{min}=4.8$, $r^{anal}_{min}=1.71$, $r_{min}=\sigma_{min}\sqrt{E_{ave}}$','Interpreter','latex','FontSize',14);
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(10055)
% subplot(2,2,1);
% plot(VM,E_TOT_vector(:,1,length(q0_vector1),1,1),'r',VM,E_TOT_vector(:,1,length(q0_vector1),1,3),'g',VM,E_TOT_vector(:,1,length(q0_vector1),1,7),'b',...
%      VM,E_TOT_vector(:,2,length(q0_vector1),1,1),'--+',VM,E_TOT_vector(:,2,length(q0_vector1),1,3),'--*',VM,E_TOT_vector(:,2,length(q0_vector1),1,7),'--o');
% xlabel('$VMs$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_{TOT}}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20, \theta=0.01','E_{ave}=20, \theta=0.25','E_{ave}=20, \theta=0.99',...
%        'E_{ave}=10, \theta=0.01','E_{ave}=10, \theta=0.25','E_{ave}=10, \theta=0.99');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $Hemogenous VMs$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(VM,E_W_vector(:,1,length(q0_vector1),1,1),'r',VM,E_W_vector(:,1,length(q0_vector1),1,3),'g',VM,E_W_vector(:,1,length(q0_vector1),1,7),'b',...
%      VM,E_W_vector(:,2,length(q0_vector1),1,1),'--+',VM,E_W_vector(:,2,length(q0_vector1),1,3),'--*',VM,E_W_vector(:,2,length(q0_vector1),1,7),'--o');
% xlabel('$VMs$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_{W}}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20, \theta=0.01','E_{ave}=20, \theta=0.25','E_{ave}=20, \theta=0.99',...
%        'E_{ave}=10, \theta=0.01','E_{ave}=10, \theta=0.25','E_{ave}=10, \theta=0.99');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $Hemogenous VMs$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(VM,r_vector(:,1,length(q0_vector1),1,1),'r',VM,r_vector(:,1,length(q0_vector1),1,3),'g',VM,r_vector(:,1,length(q0_vector1),1,7),'b',...
%      VM,r_vector(:,2,length(q0_vector1),1,1),'--+',VM,r_vector(:,2,length(q0_vector1),1,3),'--*',VM,r_vector(:,2,length(q0_vector1),1,7),'--o');
% xlabel('$VMs$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20, \theta=0.01','E_{ave}=20, \theta=0.25','E_{ave}=20, \theta=0.99',...
%        'E_{ave}=10, \theta=0.01','E_{ave}=10, \theta=0.25','E_{ave}=10, \theta=0.99');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $Hemogenous VMs$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(VM,E_TOT_vector(:,1,length(q0_vector1),1,1),'r',VM,E_TOT_vector(:,1,length(q0_vector1),1,3),'g',VM,E_TOT_vector(:,1,length(q0_vector1),1,7),'b',...
%      VM,E_TOT_vector(:,2,length(q0_vector1),1,1),'--+',VM,E_TOT_vector(:,2,length(q0_vector1),1,3),'--*',VM,E_TOT_vector(:,2,length(q0_vector1),1,7),'--o');
% xlabel('$VMs$','Interpreter','latex','FontSize',14);
% ylabel('$UTILITY$','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20, \theta=0.01','E_{ave}=20, \theta=0.25','E_{ave}=20, \theta=0.99',...
%        'E_{ave}=10, \theta=0.01','E_{ave}=10, \theta=0.25','E_{ave}=10, \theta=0.99');
% title('$q_0=N_O=240$, $s_{0}=0.01$, $Hemogenous VMs$','Interpreter','latex','FontSize',14);
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%results
% figure(100)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(:),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('without J_P, E_{ave}=0.125');
% title(['$n= $' num2str(n), ' $\theta=0.5$, $s_0= $' num2str(s(nclient,1)), ' $\sigma_{min}^{s}$= ' num2str(sigma_min_simulation),...
%        ' $\sigma_{min}^{a}$= ' num2str(sigma_min_anal),' $r_{min}^{s}$=' num2str(min(r)),...
%        ' $r_{min}^{a}$= ' num2str(r_min_anal), ' $RTT^{max}$= ' num2str(max(RTT)),...
%        ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(:),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('without J_P, E_{ave}=0.125');
% % title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
% %        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
% %        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
% %        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(:),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('without J_P, E_{ave}=0.125');
%       % title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
% %        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
% %        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
% %        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(:),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('without J_P, E_{ave}=0.125');
% % title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
% %        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
% %        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
% %        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
% grid on
% %
figure(100)
subplot(2,3,1);
plot(q0_vector1,E_W_vector(:,1),'r',q0_vector1,E_W_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
title(['$\theta=0.5$, $s_0= $' num2str(s(nclient,1)), ' $\sigma_{min}^{s}$= ' num2str(sigma_min_simulation),...
    ' $\sigma_{min}^{a}$= ' num2str(sigma_min_anal),' $r_{min}^{s}$=' num2str(min(r)),...
    ' $r_{min}^{a}$= ' num2str(r_min_anal), ' $RTT^{max}$= ' num2str(max(RTT)),...
    ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,2);
plot(q0_vector1,r_vector(:,1),'r',q0_vector1,r_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,3);
plot(q0_vector1,q_vector(:,1),'r',q0_vector1,q_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,4);
plot(q0_vector1,L_tot_vector(:,1),'r',q0_vector1,L_tot_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{L}_{tot}$ (byte)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on
%

subplot(2,3,5);
plot(q0_vector1,s_vector(:,1),'r',q0_vector1,s_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on
%

subplot(2,3,6);
plot(q0_vector1,lambda_vector(:,1),'r',q0_vector1,lambda_vector(:,2),'b');
xlabel('$q_0$','Interpreter','latex','FontSize',14);
ylabel('$\overline{\lambda}$ (byte)','Interpreter','latex','FontSize',14);
legend('J_P=0.25, E_{ave}=0.125', 'J_P=0.75, E_{ave}=0.125');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on
%


figure(10077777)
subplot(2,3,1);
plot(1:1:n,E_W(1,:),'b',1:1:n,E_W_old(1,:),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$E_{W}$ (Joule)','Interpreter','latex','FontSize',14);
legend('$E_{W}$ (Consolidate)', '$E_{W}$ (NO Consolidate)');
title(['$\overline{E}_{W}$ ' num2str(mean(E_W(1,:))), '$\overline{E}^{OLD}_{W}$ ' num2str(mean(E_W_old(1,:)))],'Interpreter','latex');
% title(['$\overline{E}_{W}$' num2str(mean(E_W(1,:))), '$J_P$=' num2str(J_P), '$E_{ave}$=' num2str(energy_ave), '$\theta=0.5$, $s_0= $' num2str(s(nclient,1)), ' $\sigma_{min}^{s}$= ' num2str(sigma_min_simulation),...
%     ' $\sigma_{min}^{a}$= ' num2str(sigma_min_anal),' $r_{min}^{s}$=' num2str(min(r)),...
%     ' $r_{min}^{a}$= ' num2str(r_min_anal), ' $RTT^{max}$= ' num2str(max(RTT)),...
%     ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on


subplot(2,3,2);
plot(1:1:n,E_TOT(:),'b',1:1:n,mean(E_TOT(:))*ones(1,n),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$E_{TOT}$ (Joule)','Interpreter','latex','FontSize',14);
legend('$E_{TOT}$', '$\overline{E}_{TOT}$');
title(['$\overline{E}_{TOT}$=' num2str(mean(E_TOT(:)))],'Interpreter','latex');
% title(['$J_P$=' num2str(J_P), '$E_{ave}$=' num2str(energy_ave), '$\theta=0.5$, $s_0= $' num2str(s(nclient,1)), ' $\sigma_{min}^{s}$= ' num2str(sigma_min_simulation),...
%     ' $\sigma_{min}^{a}$= ' num2str(sigma_min_anal),' $r_{min}^{s}$=' num2str(min(r)),...
%     ' $r_{min}^{a}$= ' num2str(r_min_anal), ' $RTT^{max}$= ' num2str(max(RTT)),...
%     ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,3);
plot(1:1:n,r(1,:),'b',1:1:n,r_old(1,:),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$r$ (byte/slot)','Interpreter','latex','FontSize',14);
legend('$r$ (Consolidate)', '$r$ (No Consolidate)');
title(['$\overline{r}$ (consolidate)= ' num2str(mean(r(1,:))), ' $\overline{r}$ (No Consolidate)= ' num2str(mean(r_old(1,:)))],'Interpreter','latex');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,4);
plot(1:1:n+1,q(1,:),'b',1:1:n+1,q_old(1,:),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$q$ (byte)','Interpreter','latex','FontSize',14);
legend('$q$ (Consolidate)', '$q$ (NO Consolidate)');
title(['$\overline{q}$ (consolidate)= ' num2str(mean(q(1,:))), ' $\overline{q}$ (No Consolidate)= ' num2str(mean(q_old(1,:)))],'Interpreter','latex');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

subplot(2,3,5);
plot(1:1:n,L_tot(1,:),'b',1:1:n,L_tot_old(1,:),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$L_{tot}$ (byte)','Interpreter','latex','FontSize',14);
legend('$L_tot$ (Consolidate)', '$L_tot$ (NO Consolidate)');
title(['$\overline{L}_{tot}$ (consolidate)= ' num2str(mean(L_tot(1,:))), ' $\overline{L}_{tot}$ (No consolidate)= ' num2str(mean(L_tot_old(1,:)))],'Interpreter','latex');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on
%

subplot(2,3,6);
plot(1:1:n+1,s(1,:),'b',1:1:n+1,s_old(1,:),'--r');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$s$ (byte)','Interpreter','latex','FontSize',14);
legend('$s$ (Consolidate)', '$s$ (NO Consolidate)');
title(['$\overline{s}$ (consolidate)= ' num2str(mean(s(1,:))), ' $\overline{s}$ (No consolidate)= ' num2str(mean(s_old(1,:)))],'Interpreter','latex');
% title(['$\theta=0.5$, $s_0=0.01$, $\sigma_{min}^{s}$=' num2str(sigma_min_simulation),...
%        '$\sigma_{min}^{a}$= ' num2str(sigma_min_anal),'$r_{min}^{s}$=' num2str(min(r)),...
%        '$r_{min}^{a}$= ' num2str(r_min_anal), '$RTT^{max}$= ' num2str(max(RTT)),...
%        '$\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on
%

figure(101)
plot(1:1:n,W_range_old(:),'-or',1:1:n,W_range(:),'b-d', 1:1:n, W_range_estimate(:),'g--s','LineWidth',2);
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$W$','Interpreter','latex','FontSize',14);
legend('$W$ (No Consolidate)', '$W$ (Consolidate)', '$W$ (Estimation+Consolidate)');
title(['$VM$= ' num2str(M), ' $ M1$= ' num2str(M1_start), '$UP_{range}$= ' num2str(up_th), '  $LOW_{range}$= ' num2str(down_th),...
    ' $E_{ave}$= ' num2str(energy_ave), ' $\theta$= ', num2str(teta),...
    ' $s_0= $' num2str(s(nclient,1)), ' $RTT^{max}$= ' num2str(max(RTT)),...
    ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

figure(102)
%subplot(2,2,4);
plot(1:1:n,iter_consolidate_slot(:),'b');
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$iteration$','Interpreter','latex','FontSize',14);
legend('$\#iter$');
title(['$\#consolidate slots$= ' num2str(sum(NO_consolidate))],'Interpreter','latex');
% title(['$J_P$=' num2str(J_P), '$E_{ave}$=' num2str(energy_ave), '$\theta=0.5$, $s_0= $' num2str(s(nclient,1)), ' $\sigma_{min}^{s}$= ' num2str(sigma_min_simulation),...
%     ' $\sigma_{min}^{a}$= ' num2str(sigma_min_anal),' $r_{min}^{s}$=' num2str(min(r)),...
%     ' $r_{min}^{a}$= ' num2str(r_min_anal), ' $RTT^{max}$= ' num2str(max(RTT)),...
%     ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

figure(103)
plot(1:1:n,active_VM_old(:),'-or',1:1:n,active_VM(:),'b-d',1:1:n, active_VM(:)+needed_new_VM(:),'g--s','LineWidth',2);
% plot(1:1:n,active_VM_old(:),'r','Marker','*','LineStyle','--',...
%      1:1:n,active_VM(:),'b','Marker','o','LineWidth',3,...
%      1:1:n, active_VM(:)+needed_new_VM(:),'g','Marker','+','LineWidth',2);
xlabel('$slot$','Interpreter','latex','FontSize',14);
ylabel('$VM actives$','Interpreter','latex','FontSize',14);
legend('Active VM (No Consolidate)', 'Active VM (consolidate)', 'Active VM (Estimate+consolidate)');
title(['$VM$= ' num2str(M), ' $ M1$= ' num2str(M1_start), '$UP_{range}$= ' num2str(up_th), '  $LOW_{range}$= ' num2str(down_th),...
    ' $E_{ave}$= ' num2str(energy_ave), ' $\theta$= ', num2str(teta),...
    ' $s_0= $' num2str(s(nclient,1)), ' $RTT^{max}$= ' num2str(max(RTT)),...
    ' $\Delta_{IP}^{max}$= ' num2str(DELTA_IP_MAX)],'Interpreter','latex');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(101)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(3,:,1),'--o',q0_vector1,E_W_vector(3,:,2),'--*', q0_vector1,E_W_vector(4,:,1),'r',q0_vector1,E_W_vector(4,:,2),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('J_P=0.25, \sigma=10','J_P=0.75, \sigma=10','J_P=0.25, \sigma=100','J_P=0.75, \sigma=100');
% title('$E_{ave}=10$ (Joule), $\theta=0.01$, $s_0=0.01$, $r_{min}=0.01$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,E_W_vector(3,:,1),'--o',q0_vector1,E_W_vector(3,:,2),'--*', q0_vector1,E_W_vector(4,:,1),'r',q0_vector1,E_W_vector(4,:,2),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('J_P=0.25, \sigma=10','J_P=0.75, \sigma=10','J_P=0.25, \sigma=100','J_P=0.75, \sigma=100');
% title('$E_{ave}=10$ (Joule), $\theta=0.01$, $s_0=0.01$, $r_{min}=0.01$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,E_W_vector(3,:,1),'--o',q0_vector1,E_W_vector(3,:,2),'--*', q0_vector1,E_W_vector(4,:,1),'r',q0_vector1,E_W_vector(4,:,2),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('J_P=0.25, \sigma=10','J_P=0.75, \sigma=10','J_P=0.25, \sigma=100','J_P=0.75, \sigma=100');
% title('$E_{ave}=10$ (Joule), $\theta=0.01$, $s_0=0.01$, $r_{min}=0.01$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,E_W_vector(3,:,1),'--o',q0_vector1,E_W_vector(3,:,2),'--*', q0_vector1,E_W_vector(4,:,1),'r',q0_vector1,E_W_vector(4,:,2),'g');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('J_P=0.25, \sigma=10','J_P=0.75, \sigma=10','J_P=0.25, \sigma=100','J_P=0.75, \sigma=100');
% title('$E_{ave}=10$ (Joule), $\theta=0.01$, $s_0=0.01$, $r_{min}=0.01$','Interpreter','latex','FontSize',14);
% grid on
% %      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(101)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(3,:,1,1),'--o',q0_vector1,E_W_vector(4,:,1,1),'--+',q0_vector1,E_W_vector(5,:,1,1),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1','E_{ave}=0.75', 'E_{ave}=0.5');
% title('$\overline{E_W}$ (Joule) , $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(3,:,1,1),'--o',q0_vector1,r_vector(4,:,1,1),'--+',q0_vector1,E_W_vector(5,:,1,1),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1','E_{ave}=0.75', 'E_{ave}=0.5');
% title('$\overline{r}$ (byte/slot), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(3,:,1,1),'--o',q0_vector1,q_vector(4,:,1,1),'--+',q0_vector1,E_W_vector(5,:,1,1),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1','E_{ave}=0.75', 'E_{ave}=0.5');
% title('$\overline{q}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(3,:,1,1),'--o',q0_vector1,s_vector(4,:,1,1),'--+',q0_vector1,E_W_vector(5,:,1,1),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1','E_{ave}=0.75', 'E_{ave}=0.5');
% title('$\overline{s}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(102)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(6,:,1,1),'--o',q0_vector1,E_W_vector(7,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{E_W}$ (Joule), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(6,:,1,1),'--o',q0_vector1,r_vector(7,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{r}$ (byte/slot), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(6,:,1,1),'--o',q0_vector1,q_vector(7,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{q}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(6,:,1,1),'--o',q0_vector1,s_vector(7,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{s}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(103)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(8,:,1,1),'--o',q0_vector1,E_W_vector(9,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{E_W}$ (Joule), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(8,:,1,1),'--o',q0_vector1,r_vector(9,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{r}$ (byte/slot), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(8,:,1,1),'--o',q0_vector1,q_vector(9,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{q}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(8,:,1,1),'--o',q0_vector1,s_vector(9,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{s}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(104)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(10,:,1,1),'--o',q0_vector1,E_W_vector(11,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{E_W}$ (Joule), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(10,:,1,1),'--o',q0_vector1,r_vector(11,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{r}$ (byte/slot), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(10,:,1,1),'--o',q0_vector1,q_vector(11,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{q}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(10,:,1,1),'--o',q0_vector1,s_vector(11,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{s}$ (byte), $\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Muinf
% figure(1000)
% subplot(2,2,1);
% plot(q0_vector1,Muinf_vector(1,:,1,1),'--o',q0_vector1,Muinf_vector(2,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=20','E_{ave}=10');
% title('$\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,Muinf_vector(3,:,1,1),'--o',q0_vector1,Muinf_vector(4,:,1,1),'--+', q0_vector1,Muinf_vector(5,:,1,1),'--o', q0_vector1,Muinf_vector(6,:,1,1),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=1', 'E_{ave}=0.75','E_{ave}=0.5', 'E_{ave}=0.25');
% title('$\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,Muinf_vector(7,:,1,1),'--o',q0_vector1,Muinf_vector(8,:,1,1),'--+', q0_vector1,Muinf_vector(9,:,1,1),'--o');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=0.125','E_{ave}=0.0650', 'E_{ave}=0.0325');
% title('$\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,Muinf_vector(10,:,1,1),'--o',q0_vector1,Muinf_vector(11,:,1,1),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\theta=0.01$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(105)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(1,:,1,2),'--o',q0_vector1,E_W_vector(2,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20','E_{ave}=10');
% title('$\overline{E_W}$ (Joule), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(1,:,1,2),'--o',q0_vector1,r_vector(2,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20','E_{ave}=10');
% title('$\overline{r}$ (byte/slot), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(1,:,1,2),'--o',q0_vector1,q_vector(2,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20','E_{ave}=10');
% title('$\overline{q}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(1,:,1,2),'--o',q0_vector1,s_vector(2,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=20','E_{ave}=10');
% title('$\overline{s}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% %      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(106)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(3,:,1,2),'--o',q0_vector1,E_W_vector(4,:,1,2),'--+',q0_vector1,E_W_vector(5,:,1,2),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1', 'E_{ave}=0.75','E_{ave}=0.5');
% title('$\overline{E_W}$ (Joule), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(3,:,1,2),'--o',q0_vector1,r_vector(4,:,1,2),'--+',q0_vector1,E_W_vector(5,:,1,2),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1', 'E_{ave}=0.75','E_{ave}=0.5');
% title('$\overline{r}$ (byte/slot), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(3,:,1,2),'--o',q0_vector1,q_vector(4,:,1,2),'--+',q0_vector1,E_W_vector(5,:,1,2),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1', 'E_{ave}=0.75','E_{ave}=0.5');
% title('$\overline{q}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(3,:,1,2),'--o',q0_vector1,s_vector(4,:,1,2),'--+',q0_vector1,E_W_vector(5,:,1,2),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=1', 'E_{ave}=0.75','E_{ave}=0.5');
% title('$\overline{s}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(107)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(6,:,1,2),'--o',q0_vector1,E_W_vector(7,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{E_W}$ (Joule), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(6,:,1,2),'--o',q0_vector1,r_vector(7,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{r}$ (byte/slot), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(6,:,1,2),'--o',q0_vector1,q_vector(7,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{q}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(6,:,1,2),'--o',q0_vector1,s_vector(7,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.25','E_{ave}=0.125');
% title('$\overline{s}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(108)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(8,:,1,2),'--o',q0_vector1,E_W_vector(9,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{E_W}$ (Joule), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(8,:,1,2),'--o',q0_vector1,r_vector(9,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{r}$ (byte/slot), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(8,:,1,2),'--o',q0_vector1,q_vector(9,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{q}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(8,:,1,2),'--o',q0_vector1,s_vector(9,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.0650','E_{ave}=0.0325');
% title('$\overline{s}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(109)
% subplot(2,2,1);
% plot(q0_vector1,E_W_vector(10,:,1,2),'--o',q0_vector1,E_W_vector(11,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{E_W}$ (Joule)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{E_W}$ (Joule), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,r_vector(10,:,1,2),'--o',q0_vector1,r_vector(11,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{r}$ (byte/slot)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{r}$ (byte/slot), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,q_vector(10,:,1,2),'--o',q0_vector1,q_vector(11,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{q}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{q}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,s_vector(10,:,1,2),'--o',q0_vector1,s_vector(11,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('$\overline{s}$ (byte)','Interpreter','latex','FontSize',14);
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\overline{s}$ (byte), $\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Muinf
% figure(1001)
% subplot(2,2,1);
% plot(q0_vector1,Muinf_vector(1,:,1,2),'--o',q0_vector1,Muinf_vector(2,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=20','E_{ave}=10');
% title('$\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,2);
% plot(q0_vector1,Muinf_vector(3,:,1,2),'--o',q0_vector1,Muinf_vector(4,:,1,2),'--+', q0_vector1,Muinf_vector(5,:,1,2),'--o', q0_vector1,Muinf_vector(6,:,1,2),'--*');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=1','E_{ave}=0.75','E_{ave}=0.5', 'E_{ave}=0.25');
% title('$\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,3);
% plot(q0_vector1,Muinf_vector(7,:,1,2),'--o',q0_vector1,Muinf_vector(8,:,1,2),'--+', q0_vector1,Muinf_vector(9,:,1,2),'--o');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=0.125','E_{ave}=0.0650', 'E_{ave}=0.0325');
% title('$\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on
%
% subplot(2,2,4);
% plot(q0_vector1,Muinf_vector(10,:,1,2),'--o',q0_vector1,Muinf_vector(11,:,1,2),'--+');
% xlabel('$q_0$','Interpreter','latex','FontSize',14);
% ylabel('Mu_{inf}');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% title('$\theta=0.99$, $s_0=240$','Interpreter','latex','FontSize',14);
% grid on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lambda
% figure(1002)
% subplot(2,2,1);
% plot([1, 5, 8, 10, 25, 50, 100],lambda_vector(1,:),'--o',[1, 5, 8, 10, 25, 50, 100],lambda_vector(2,:),'--+');
% xlabel('q_0');
% ylabel('mean(\lambda)');
% legend('E_{ave}=1','E_{ave}=0.5');
% grid on
%
% subplot(2,2,2);
% plot([1, 5, 8, 10, 25, 50, 100],lambda_vector(3,:),'--o',[1, 5, 8, 10, 25, 50, 100],lambda_vector(4,:),'--+');
% xlabel('q_0');
% ylabel('mean(\lambda)');
% legend('E_{ave}=0.25','E_{ave}=0.125');
% grid on
%
% subplot(2,2,3);
% plot([1, 5, 8, 10, 25, 50, 100],lambda_vector(5,:),'--o',[1, 5, 8, 10, 25, 50, 100],lambda_vector(6,:),'--+');
% xlabel('q_0');
% ylabel('mean(\lambda)');
% legend('E_{ave}=0.125/2','E_{ave}=0.125/4');
% grid on
%
% subplot(2,2,4);
% plot([1, 5, 8, 10, 25, 50, 100],lambda_vector(7,:),'--o',[1, 5, 8, 10, 25, 50, 100],lambda_vector(8,:),'--+');
% xlabel('q_0');
% ylabel('mean(\lambda)');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% r
% figure(1003)
% subplot(2,2,1);
% plot([1, 5, 8, 10, 25, 50, 100],r_vector(1,:),'--o',[1, 5, 8, 10, 25, 50, 100],r_vector(2,:),'--+');
% xlabel('q_0');
% ylabel('mean(r^*)');
% legend('E_{ave}=1','E_{ave}=0.5');
% grid on
%
% subplot(2,2,2);
% plot([1, 5, 8, 10, 25, 50, 100],r_vector(3,:),'--o',[1, 5, 8, 10, 25, 50, 100],r_vector(4,:),'--+');
% xlabel('q_0');
% ylabel('mean(r^*)');
% legend('E_{ave}=0.25','E_{ave}=0.125');
% grid on
%
% subplot(2,2,3);
% plot([1, 5, 8, 10, 25, 50, 100],r_vector(5,:),'--o',[1, 5, 8, 10, 25, 50, 100],r_vector(6,:),'--+');
% xlabel('q_0');
% ylabel('mean(r^*)');
% legend('E_{ave}=0.125/2','E_{ave}=0.125/4');
% grid on
%
% subplot(2,2,4);
% plot([1, 5, 8, 10, 25, 50, 100],r_vector(7,:),'--o',[1, 5, 8, 10, 25, 50, 100],r_vector(8,:),'--+');
% xlabel('q_0');
% ylabel('mean(r^*)');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L_tot
% figure(1004)
% subplot(2,2,1);
% plot([1, 5, 8, 10, 25, 50, 100],L_tot_vector(1,:),'--o',[1, 5, 8, 10, 25, 50, 100],L_tot_vector(2,:),'--+');
% xlabel('q_0');
% ylabel('mean(L_{tot}) (byte)');
% legend('E_{ave}=1','E_{ave}=0.5');
% grid on
%
% subplot(2,2,2);
% plot([1, 5, 8, 10, 25, 50, 100],L_tot_vector(3,:),'--o',[1, 5, 8, 10, 25, 50, 100],L_tot_vector(4,:),'--+');
% xlabel('q_0');
% ylabel('mean(L_{tot}) (byte)');
% legend('E_{ave}=0.25','E_{ave}=0.125');
% grid on
%
% subplot(2,2,3);
% plot([1, 5, 8, 10, 25, 50, 100],L_tot_vector(5,:),'--o',[1, 5, 8, 10, 25, 50, 100],L_tot_vector(6,:),'--+');
% xlabel('q_0');
% ylabel('mean(L_{tot}) (byte)');
% legend('E_{ave}=0.125/2','E_{ave}=0.125/4');
% grid on
%
% subplot(2,2,4);
% plot([1, 5, 8, 10, 25, 50, 100],r_vector(7,:),'--o',[1, 5, 8, 10, 25, 50, 100],r_vector(8,:),'--+');
% xlabel('q_0');
% ylabel('mean(L_{tot}) (byte)');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% s
% figure(1005)
% subplot(2,2,1);
% plot([1, 5, 8, 10, 25, 50, 100],s_vector(1,:),'--o',[1, 5, 8, 10, 25, 50, 100],s_vector(2,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{s}$ (byte)');
% legend('E_{ave}=1','E_{ave}=0.5');
% grid on
%
% subplot(2,2,2);
% plot([1, 5, 8, 10, 25, 50, 100],s_vector(3,:),'--o',[1, 5, 8, 10, 25, 50, 100],s_vector(4,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{s}$ (byte)');
% legend('E_{ave}=0.25','E_{ave}=0.125');
% grid on
%
% subplot(2,2,3);
% plot([1, 5, 8, 10, 25, 50, 100],s_vector(5,:),'--o',[1, 5, 8, 10, 25, 50, 100],s_vector(6,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{s}$ (byte)');
% legend('E_{ave}=0.125/2','E_{ave}=0.125/4');
% grid on
%
% subplot(2,2,4);
% plot([1, 5, 8, 10, 25, 50, 100],s_vector(7,:),'--o',[1, 5, 8, 10, 25, 50, 100],s_vector(8,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{s}$ (byte)');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% grid on
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% q
% figure(1006)
% subplot(2,2,1);
% plot([1, 5, 8, 10, 25, 50, 100],q_vector(1,:),'--o',[1, 5, 8, 10, 25, 50, 100],q_vector(2,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{q}$ (byte)');
% legend('E_{ave}=1','E_{ave}=0.5');
% grid on
%
% subplot(2,2,2);
% plot([1, 5, 8, 10, 25, 50, 100],q_vector(3,:),'--o',[1, 5, 8, 10, 25, 50, 100],q_vector(4,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{q}$ (byte)');
% legend('E_{ave}=0.25','E_{ave}=0.125');
% grid on
%
% subplot(2,2,3);
% plot([1, 5, 8, 10, 25, 50, 100],q_vector(5,:),'--o',[1, 5, 8, 10, 25, 50, 100],q_vector(6,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{q}$ (byte)');
% legend('E_{ave}=0.125/2','E_{ave}=0.125/4');
% grid on
%
% subplot(2,2,4);
% plot([1, 5, 8, 10, 25, 50, 100],q_vector(7,:),'--o',[1, 5, 8, 10, 25, 50, 100],q_vector(8,:),'--+');
% xlabel('q_0');
% ylabel('$\overline{q}$ (byte)');
% legend('E_{ave}=0.01','E_{ave}=0.005');
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%real results
%  figure(20)
% figure(1000)
% plot(1:1:n,E_W(1,:),'--o',1:1:n,mean(E_W(1,:)).*ones(1,n),'--+');
% xlabel('slots');
% ylabel('E_W (Joule)');
% legend('E_W','$\overline{E_W}$');
%
%
% figure(1001)
% plot(1:1:n+1,Mu(1,:),'--o',1:1:n+1,mean(Mu(1,:)).*ones(1,n+1),'--+');
% xlabel('slots');
% ylabel('Mu');
% legend('Mu','mean(Mu)');
%
%
% figure(1002)
% plot(1:1:n,r(1,:),'--o',1:1:n,mean(r(1,:)).*ones(1,n),'--+');
% xlabel('slots');
% ylabel('r (byte)');
% legend('r','$\overline{r}$');
%
% figure(1003)
% plot(1:1:n,lambda(1,:),'--o',1:1:n,mean(lambda(1,:)).*ones(1,n),'--+');
% xlabel('slots');
% ylabel('lambda (byte)');
% legend('lambda','mean(lambda)');
%
%
% figure(1004)
% plot(1:1:n+1,q(1,:),'--o',1:1:n+1,mean(q(1,:)).*ones(1,n+1),'--+');
% xlabel('slots');
% ylabel('q (byte)');
% legend('q','$\overline{q}$');
%
%
% figure(1005)
% plot(1:1:n+1,s(1,:),'--o',1:1:n+1,mean(s(1,:)).*ones(1,n+1),'--+');
% xlabel('slots');
% ylabel('s (byte)');
% legend('s','$\overline{s}$');
%
% figure(1006)
% plot(1:1:n,sigma(1,:),'--o',1:1:n,mean(sigma(1,:)).*ones(1,n),'--+');
% xlabel('slots');
% ylabel('sigma (byte)');
% legend('\sigma','mean(\sigma)');
%
% figure(1007)
% plot(1:1:n,L_tot(:),'--o',1:1:n,mean(L_tot(:)).*ones(1,n),'--+');
% xlabel('slots');
% ylabel('L_{tot} (byte)');
% legend('L_{tot}','mean(L_{tot})');
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DVFS testing
% figure(3000+NumSimulazioni+1)
% plot(VM,costo_tot_VM,'--o',VM,costo_tot_VM_f_succ,'--o',VM,costo_tot_VM_f_prec,'--o',VM,costo_tot_VM_TS,'--o',VM,costo_tot_VM_TSpiuSwitch,'--o'),xlabel('VM'),ylabel('Overall Cost'),legend('Q=\infty (no DVFS)','f_succ','f_prec','TS','TS+Switch');
%
% figure(3000+NumSimulazioni+2)
% plot(VM,costo_VM_TS_switch,'--o'),xlabel('VM'),ylabel('Time Sharing contiguous switch cost');
%
% figure(3000+NumSimulazioni+3)
% plot(VM,f_VM_prima,'--o',VM,f_VM_ultima,'--o'),xlabel('VM'),ylabel('rate VM(1)'),legend('first VM','last VM');
%
% %figure(21)
% figure(3000+NumSimulazioni+4)
% plot(VM,L_VM_prima,'--o',VM,L_VM_ultima,'--o'),xlabel('VM'),ylabel('workload VM(1)'),legend('first VM','last VM');









