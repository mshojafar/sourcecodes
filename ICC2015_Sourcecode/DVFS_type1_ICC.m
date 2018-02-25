%ScriptMax 1a cond + ni cond blocco totale version 4
close all;
clear all;
clc;

tic;
VMNo=100;                                        % VM numbers
N = 2000;                                        %number of workloads

nVM=[];
% +-+-+-+-+-+ parameters +-+-+-+-+-+ 
for i=1:VMNo
nVM = [nVM,i];                        %number of virtual machines
end


speed_cost_all = zeros(length(nVM),N);           %matrix of speed cost_all (1 x workloads)
switch_cost_all = zeros(length(nVM),N);          %matrix of switch cost (1 x workloads)
channel_cost_all = zeros(length(nVM),N);         %matrix of channel cost (1 x workloads)
tot_cost_all = zeros(length(nVM),N);             %matrix of total cost (1 x workloads)


k_e = 0.005;                            %frequency change cost  parameter
delta = 0.1;                          %task time
f0 = 0;                             %percentage of fmax of initial frequency
iterations = 10000;                 %max number of iteration
gamma = 0.5;

tol_workload_allocated = 0.01;      %tolerance 1st condition
epsilon = 0.01;                     %tolerance 2nd condition
R_tot = 100;                        %capacity of NetDC
%L_tot_vect = 6+ round((10-6) * rand(1,N)); % PMR=1.25, L_tot=8+-2
%save('workloads8PMR125No1000.mat', 'L_tot_vect');                  %loading 100 fixed random workloads
L_tot_vect=load('workloads8PMR125No2000.mat','L_tot_vect');                                %loading 100 fixed random workloads
L_tot_vect = L_tot_vect.L_tot_vect;

T_tot = 5;      
beta = 0.41/(T_tot - delta)^1.995;%time out
%beta = 0.34/(T_tot - delta)^1;%time out

rule_counter=zeros(1,N);                    % (insit of notsatisfactin cases) states for the each incoming workload facing in 3 rules
rule1_counter=0;
rule2_counter=0;
rule3_counter=0;

Zeta_vector=[0.5:0.1:100];
%VMs=[2, 100];

% iter_mat1 = zeros(1,N); 
% mu_mat1 = zeros(N,iterations);
% ni_mat1 = zeros(N,iterations);
% f_opt_mat1=zeros(N,iterations);
% L_mat1=zeros(N,iterations);

%for M = nVM(1):VMNo
for M = 100:100
%for M = 1:length(VMs)
    %M=VMs(M);
    M
    f_max = 105 * ones(1,M);        % (Mbit/s) max working frequency
    f_zero = (f0/100) * f_max;      %initial working frequency
    %f_zero = 0;      %initial working frequency
    W = ones(1,M);                  %band of each channel (MHz)
    %Zeta = 0.5 * ones(1,M);         %channel parameters (mW)
    Zeta = Zeta_vector(1:M); 
    f_opt_new=zeros(1,M);
%     N_0 = ones(1,M);g = ones(1,M);Zeta = (N_0 .* W) / g;R = [0:0.1:R_max];P_net = Zeta .* (2^(R/W) - 1);
    E_max = 60 * ones(1,M);         %max power consumption (@fmax)
    omega = ones(1,M);              %relative power consumption (virtualization cost)
    Delta = delta *  ones(1,M);     %task time
    Delta_max = max(Delta);         %max task time
    L_b = zeros(1,M);               %background task size for each VM
    
%    speed_cost = zeros(N,iterations);           %matrix of speed cost (workloads x iteration)
%    switch_cost = zeros(N,iterations);          %matrix of switch cost (workloads x iteration)
%    channel_cost = zeros(N,iterations);         %matrix of channel cost (workloads x iteration)
%    tot_cost = zeros(N,iterations);             %matrix of total cost (workloads x iteration)
    
        
    temp=2.*k_e+(2.*E_max.*omega./(f_max.^2));  %temp var for calculation of fi*
    Th = 2 * (Zeta./W) * log(2);                %threshold = derivation of Pnet(Ri) 
    
    iter_mat = zeros(1,N);                      %number of iteration necessary for 1st condition for each workload
    
%     MatrixRis = zeros(1+M+M,N);                 %matrix of results, row: 1=final_iter 2=final_alpha 3=tot_cost 4=final_ni 5=2nd_condition 6=(ni<esp & g(.)<esp)
%     L_allocated = zeros(N,iterations);          %matrix of total allocation for each iteration fon each workload 
%     f_opt_mat = zeros(N,iterations);            %matrix of frequency allocation for each iteretion for each workload FOR THE 1ST VM
%     L_mat = zeros(N,iterations);                %matrix of  allocation for each iteretion for each workload FOR THE 1ST VM
%     ni_mat = zeros(N,iterations);               %matrix of "ni" value for each iteretion for each workload FOR THE 1ST VM
%     mu_mat = zeros(N,iterations);
%     
    %matrix of "mu" value for each iteretion for each workload
    alpha_mat = zeros(N,iterations);            %matrix of "alpha" value for each iteretion for each workload
    pippo_mat = zeros(N,iterations);            %matrix of "pippo" value for each iteretion for each workload FOR THE 1ST VM
    
    for l = 1:N                                 % workload number
        %check feasibility
        feas_workload = (sum(f_max.*Delta-L_b) >= L_tot_vect(l));
        feas_time = (L_tot_vect(l) <= R_tot.*(T_tot-Delta_max)./2);
        feas_back = min(Delta.*f_max >= L_b);
        feasibility = feas_workload && feas_time && feas_back;
        if ~ feasibility
            error('Problem unfeasible for the workload %d',l);
        else
            fprintf('Problem feasibile for the workload %d\n',l);
            Delta_load = zeros(1,iterations);
            L = zeros(iterations,M);
            alpha = zeros(1,iterations);
            pippo = zeros(1,iterations);
            V = zeros(1,iterations);
            V_pippo = zeros(1,iterations);
            mu = zeros(1,iterations);
            ni = zeros(iterations,M);
            y = zeros(iterations,M);
            f_opt = zeros(iterations,M);
            i = 2;
            while i <= iterations
                bandierina = 0;
                k=0;
                cond3 = zeros(1,M);
                cond2 = zeros(1,M);
                %calculate of MU
                mu(i) = max(0,(mu(i-1) - ((alpha(i-1))) * (sum(L(i-1,:)) - L_tot_vect(l)) ) );
                
                %calculate of temp variable
                y(i,:) = max(0,( ((T_tot-Delta)/2) .* W .* log2(mu(i) ./ Th) ) );
                
                %calculate of ni; the old was: ni(i,:) = max(0,( ni(i-1,:) + (alpha(i-1) * ( L(i-1,:) - (Delta.*f_opt(i-1,:))))))
                ni(i,:) = (ni(i-1,:) + (pippo(i-1) * (y(i,:) - Delta.*f_opt(i-1,:))));

                %calculate of f_opt
                f_star = (2*k_e*f_zero + ni(i,:).*Delta)./ temp;
                f_opt(i,:) = max(L_b/Delta,min(f_star,f_max)); 
                               
                %calculate of L
                L(i,:) = min(f_opt(i,:) .* Delta(1,:),y(i,:));

%                 ni(i,ni(i,:)<=0.01) = 0;
%                 for h=1:M
%                     if ni(i,h)>0
%                         L(i,h)=f_opt(i,h) .* Delta(1,h);
%                     else
%                         L(i,h) = y(i,h);
%                     end
%                 end
                
                % 2nd and 3rd condition
                for x = 1:M
                    cond3(1,x) = (((abs(ni(i,x)))<=epsilon)||((abs(L(i,x)-f_opt(i,x)*Delta(1,x)))<=epsilon));
                    cond2(1,x) = L(i,x)<=f_opt(i,x)*Delta(1,x);
                end
                
%                 if (  (sum(cond3(1,:)) < M)||(sum(cond2(1,:)) < M)  )
%                     alpha(i-1) = 0.5 * alpha(i-1);
%                     bandierina = 1;
%                     k = k + 1;
%                     fprintf('alpha loop...);
%                 end  
                if (bandierina == 0 && k<100)
                    alpha(i) = max(0,min(beta,alpha(i-1) - gamma * V(i-1) * (sum(L(i-1,:)) - L_tot_vect(l))));
                    V(i) = (1 - alpha(i-1)) * V(i-1) - (sum(L(i-1,:)) - L_tot_vect(l));
                    pippo(i) = max(0,min(0.001,pippo(i-1) - gamma * V_pippo(i-1) *(Delta(1)*f_opt(i-1,1)- y(i,1))));
                    V_pippo(i) = (1 - pippo(i-1)) * V_pippo(i-1) - (Delta(1)*f_opt(i-1,1)- y(i,1));
                    Delta_load(i)=sum(L(i,:))-L_tot_vect(l);
%                     speed_cost(l,i) = sum(((f_opt(i,:)./f_max).^2).*omega.*E_max);                                % Computing (11.1)
%                     switch_cost(l,i) = sum(k_e.*(f_opt(i,:)-f_zero).^2);                                          % Switching (11.1)
%                     channel_cost(l,i) = (T_tot-Delta(1)).*sum(Zeta .* (2.^(2*L(i,:)./((T_tot-Delta).*W))-1));     % Networking
%                     tot_cost(l,i) = speed_cost(l,i) + switch_cost(l,i) + channel_cost(l,i);                       % Overall
%                     
%                    L_allocated(l,i) = sum(L(i,:));
%                     L_mat(l,i) = L(i,1);
%                     f_opt_mat(l,i) = f_opt(i,1);
%                     ni_mat(l,i) = ni(i,1);
%                     mu_mat(l,i) = mu(i);
%                     alpha_mat(l,i) = alpha(i);
%                     pippo_mat(l,i) = pippo(i);
                     cond_car_all = abs(Delta_load(i)/L_tot_vect(l))>=tol_workload_allocated;
                    if cond_car_all
                        iter_mat(l) = i;
%                         MatrixRis(1:M,l) = L(i,:);
%                         MatrixRis(M+1:M+M,l) = f_opt(i,:);
%                         MatrixRis(M+M+1,l) = tot_cost(l,i);
                    else
                        fprintf('correct allocation for workload n: %d (%d) at iteration %d\n',l,L_tot_vect(l),i);
                        iter_mat(l) = i;
%                         MatrixRis(1:M,l) = L(i,:);
%                         MatrixRis(M+1:M+M,l) = f_opt(i,:);
%                         MatrixRis(M+M+1,l) = tot_cost(l,i);
                        break;
                    end
                    i = i + 1;
                end
            end
        end
        if i == iterations + 1
            i = i - 1;
        end
        
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DESCRETE
        
        f_precedente=zeros(1,M);
        f_successiva=zeros(1,M);
        x1=zeros(1,M);
        
        Q=6; %numero frequenze discrete;
        f_discrete=zeros(1,Q);
        for conta=1:M
            %f_discrete(conta,:)=[0:f_max(conta)./(Q-1):f_max(conta)];
            f_discrete(conta,:)=[0, 5,50,70,90,105]; % Mhz or Mb/s
        end
        
        f_opt_new=f_opt(iter_mat(l), :);
        
        for conta=1:M
            
            delta_f_discrete=f_discrete(conta,:)-f_opt_new(conta);
            [ff,ind_ff]=min(abs(delta_f_discrete));
            if ff==0
                f_precedente(conta)=f_discrete(conta,ind_ff);
                f_successiva(conta)=f_discrete(conta,ind_ff);
                x1(conta)=1; %qualsiasi valore è indifferente
            elseif ind_ff==1
                f_precedente(conta)=f_discrete(conta,1);
                f_successiva(conta)=f_discrete(conta,2);
                x1(conta)=abs(f_opt_new(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
            elseif ind_ff==Q
                f_precedente(conta)=f_discrete(conta,Q-1);
                f_successiva(conta)=f_discrete(conta,Q);
                x1(conta)=abs(f_opt_new(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
            elseif delta_f_discrete(ind_ff)>0
                f_precedente(conta)=f_discrete(conta,ind_ff-1);
                f_successiva(conta)=f_discrete(conta,ind_ff);
                x1(conta)=abs(f_opt_new(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
            else
                f_precedente(conta)=f_discrete(conta,ind_ff);
                f_successiva(conta)=f_discrete(conta,ind_ff+1);
                x1(conta)=abs(f_opt_new(conta)-f_precedente(conta))./(f_successiva(conta)-f_precedente(conta));
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         %speed_cost_all(M,l)=speed_cost(l,iter_mat(l));
       % switch_cost_all(M,l)=switch_cost(l,iter_mat(l));
        
%         speed_cost_all(M,l)=sum(((((f_precedente/f_max).^2)).*omega.*E_max).*(1-x1)+((((f_successiva/f_max).^2)).*omega.*E_max).*x1);
%         switch_cost_all(M,l)=sum((k_e.*(f_precedente-f_zero).^2).*(1-x1)+(k_e.*(f_successiva-f_zero).^2).*x1);
%         
%         channel_cost_all(M,l)=channel_cost(l,iter_mat(l));
%         
%         tot_cost_all(M,l)=speed_cost_all(M,l)+switch_cost_all(M,l)+channel_cost_all(M,l);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        f_zero = f_opt(i,:);
        %1a Condizione
        if i > iterations
            i=iterations;
        end
        if i >= iterations
            fprintf('!!!!!!!WARNING!!!!!!!prima condizione non soddifatta: NON CONVERGE!!!\n');
            rule_counter(l)=1;
            rule1_counter=rule1_counter+1;
        end
        %3a Condizione
        for x = 1:M
            cond3(1,x) = (((abs(ni(i-1,x)))<=epsilon)||((abs(L(i-1,x)-f_opt(i-1,x)*Delta(1,x)))<=epsilon));
        end
        if sum(cond3(1,:))< M
            fprintf('!!!!!!!WARNING!!!!!!!terza condizione NON soddisfatta\n');
            rule_counter(l)=3;
            rule3_counter=rule3_counter+1;
        end
        %2a condizione
        if L(i,:)>(f_opt(i,:) .* Delta(1,:));
            fprintf('!!!!!!!WARNING!!!!!!!seconda condizione NON soddisfatta\n'); 
            rule_counter(l)=2;
            rule2_counter=rule2_counter+1;            
        end
    end
    clear speed_cost switch_cost channel_cost tot_cost 
%     if M==2
%     mu_mat1 = mu_mat;
%     iter_mat1=iter_mat;
%     ni_mat1 = ni_mat;
%     f_opt_mat1=f_opt_mat;
%     L_mat1=L_mat;
%     end
end
time = toc;
fprintf('tempo: %d\n',time);

NNN=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DESCRETE COMPARISONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Presimulation-definitions
% speed_cost_all_WL=zeros(1,N);
% switch_cost_all_WL=zeros(1,N);
% channel_cost_all_WL=zeros(1,N);
% tot_cost_all_WL=zeros(1,N);
% 
% speed_cost_all_VM=zeros(1,length(nVM));
% switch_cost_all_VM=zeros(1,length(nVM));
% channel_cost_all_VM=zeros(1,length(nVM));
% tot_cost_all_VM=zeros(1,length(nVM));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Presimulation-for
% for j=1:N                   % workload
%     speed_cost_all_WL(j)=mean(speed_cost_all(:,j));
%     switch_cost_all_WL(j)=mean(switch_cost_all(:,j));
%     channel_cost_all_WL(j)=mean(channel_cost_all(:,j));
%     tot_cost_all_WL(j)=mean(tot_cost_all(:,j));
% end
% 
% for j=1:length(nVM)             % VMs 
%     speed_cost_all_VM(j)=mean(speed_cost_all(j,:));
%     switch_cost_all_VM(j)=mean(switch_cost_all(j,:));
%     channel_cost_all_VM(j)=mean(channel_cost_all(j,:));
%     tot_cost_all_VM(j)=mean(tot_cost_all(j,:));
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % figure(100)                     % non staisfaction numbers
% % plot(,rule1_counter,'--+');
% % xlabel('Workload');
% % ylabel('$\overline{\mathcal{E}}_{CPU}$');
% % grid on
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(100)                     % speed
% plot(1:1:N,speed_cost_all_WL(:));
% xlabel('Workload');
% ylabel('$\overline{\mathcal{E}}_{CPU}$');
% grid on
% 
% figure(101)                     % switch
% plot(1:1:N,switch_cost_all_WL(:));
% xlabel('Workload');
% ylabel('$\overline{\mathcal{E}}_{switch}$');
% grid on
% 
% figure(102)                     % channel
% plot(1:1:N,channel_cost_all_WL(:));
% xlabel('Workload');
% ylabel('$\overline{\mathcal{E}}^{net}$');
% grid on
% 
% figure(103)                     % total
% plot(1:1:N,tot_cost_all_WL(:));
% xlabel('Workload');
% ylabel('$\overline{\mathcal{E}}^{tot}$');
% grid on
% 
% figure(104)                     % speed
% plot(1:1:VMNo,speed_cost_all_VM(:));
% xlabel('VM');
% ylabel('$\overline{\mathcal{E}}_{CPU}$');
% grid on
% 
% figure(105)                     % switch
% plot(1:1:VMNo,switch_cost_all_VM(:));
% xlabel('VM');
% ylabel('$\overline{\mathcal{E}}_{switch}$');
% grid on
% 
% figure(106)                     % channel
% plot(1:1:VMNo,channel_cost_all_VM(:));
% xlabel('VM');
% ylabel('$\overline{\mathcal{E}}^{net}$');
% grid on
% 
% figure(107)                     % total
% plot(1:1:VMNo,tot_cost_all_VM(:));
% xlabel('VM');
% ylabel('$\overline{\mathcal{E}}^{tot}$');
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu_last=zeros (1,N);
% ni_last=zeros (1,N);
% alpha_last=zeros (1,N);
% pippo_last=zeros (1,N);
% 
% for i=1:N
% mu_last(i)=mu_mat(i,iter_mat(i));
% ni_last(i)=ni_mat(i,iter_mat(i));
% alpha_last(i)=alpha_mat(i,iter_mat(i));
% pippo_last(i)=pippo_mat(i,iter_mat(i));
% end
% 
% figure(1000)
% plot(1:1:N,mu_last(:),'r');
% title('mu');
% xlabel('WL');
% ylabel('$\mu$');
% grid on
% 
% figure(1001)
% plot(1:1:N,ni_last(:),'r');
% title('ni');
% xlabel('WL');
% ylabel('$\ni$');
% grid on
% 
% figure(1022)
% plot(1:1:iter_mat1(50),mu_mat1(50,1:iter_mat1(50)),'r',1:1:iter_mat1(1000),mu_mat1(1000,1:iter_mat1(1000)),'b',1:1:iter_mat(50),mu_mat(50,1:iter_mat(50)),'r',1:1:iter_mat(1000),mu_mat(1000,1:iter_mat(1000)),'b');
% title('$\mu$ in $\Delta=0.1$ (s), $\gamma=0.5$, $T_t=5$ (s), $\beta=0.15$');
% legend('M=2, $\delta=0.1$, N=50','M=2, $\delta=0.1$, N=1000', 'M=100, $\delta=0.1$, N=50','M=100, $\delta=0.1$, N=1000')
% xlabel('iterations');
% ylabel('$\mu$');
% grid on
% 
% figure(1023)
% plot(1:1:iter_mat1(50),ni_mat1(50,1:iter_mat1(50)),'r',1:1:iter_mat1(1000),ni_mat1(1000,1:iter_mat1(1000)),'b',1:1:iter_mat(50),ni_mat(50,1:iter_mat(50)),'r',1:1:iter_mat(1000),ni_mat(1000,1:iter_mat(1000)),'b');
% title('$\ni$ in $\Delta=0.1$ (s), $\gamma=0.5$, $T_t=5$ (s), $\beta=0.15$');
% legend('M=2, N=50','M=2, N=1000', 'M=100, N=50','M=100, N=1000')
% xlabel('iterations');
% ylabel('$\ni$');
% grid on
% 
% figure(1024)
% plot(1:1:iter_mat1(50),f_opt_mat1(50,1:iter_mat1(50)),'r',1:1:iter_mat1(1000),f_opt_mat1(1000,1:iter_mat1(1000)),'b',1:1:iter_mat(50),f_opt_mat(50,1:iter_mat(50)),'r',1:1:iter_mat(1000),f_opt_mat(1000,1:iter_mat(1000)),'b');
% title('$f_{opt}$ in $\Delta=0.1$ (s), $k_e=0.005$, $f_i^{max}=105$, $\gamma=0.5$, $T_t=5$ (s), $\beta=0.15$, $\overline{\mathcal{E}}_i^{max}=60$ $(J)$');
% legend('M=2, N=50','M=2, N=1000', 'M=100, N=50','M=100, N=1000')
% xlabel('iterations');
% ylabel('$f_{opt}$');
% grid on
% 
% figure(1025)
% plot(1:1:iter_mat1(50),L_mat1(50,1:iter_mat1(50)),'r',1:1:iter_mat1(1000),L_mat1(1000,1:iter_mat1(1000)),'b',1:1:iter_mat(50),L_mat(50,1:iter_mat(50)),'r',1:1:iter_mat(1000),L_mat(1000,1:iter_mat(1000)),'b');
% title('$L_{opt}$ in $\Delta=0.1$ (s), $k_e=0.005$, $f_i^{max}=105$, $\gamma=0.5$, $T_t=5$ (s), $\beta=0.15$, $\overline{\mathcal{E}}_i^{max}=60$ $(J)$');
% legend('M=2, N=50','M=2, N=1000', 'M=100, N=50','M=100, N=1000')
% xlabel('iterations');
% ylabel('$L_{opt}$');
% grid on

% figure(1002)
% plot(1:1:N,alpha_last(:),'r');
% title('\alpha');
% xlabel('WL');
% ylabel('$\alpha$');
% grid on
% 
% figure(1003)
% plot(1:1:N,pippo_last(:),'r');
% title('\alpha_{NEW}');
% xlabel('WL');
% ylabel('$\alpha_{NEW}$');
% grid on
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure()
% plot(mu_mat(NNN,2:iter_mat(NNN)),'r');
% title('mu');
% 
% figure()
% plot(ni_mat(NNN,2:iter_mat(NNN)),'r');
% title('ni');

% figure()
% plot(L_mat(NNN,2:iter_mat(NNN)),'r');
% title('L');
% 
% figure()
% plot(f_opt_mat(NNN,2:iter_mat(NNN)),'r');
% title('f ottimo');

% figure()
% plot(alpha_mat(NNN,2:iter_mat(NNN)),'r');
% title('alpha');
% 
% figure()
% plot(pippo_mat(NNN,2:iter_mat(NNN)),'r');
% title('pippo');

% figure()
% bar(MatrixRis(1:M,NNN));
% title('allocazione carico (L) per il wl 1');
% 
% figure()
% bar(MatrixRis(1+M:M+M,NNN));
% title('allocazione frequenze (f_opt) per il wl 1');


