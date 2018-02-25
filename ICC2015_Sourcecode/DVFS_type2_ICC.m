clc
clear all
close all

VM=[1:1:100];%numero Macchine Virtuali
NumSimulazioni=length(VM);

%IDEAL
costo_speed_VM=zeros(1,NumSimulazioni);
costo_switch_VM=zeros(1,NumSimulazioni);
costo_channel_VM=zeros(1,NumSimulazioni);
costo_tot_VM=zeros(1,NumSimulazioni);
%NICOLA (COMNET)
speed_cost_all_VM=zeros(1,NumSimulazioni);
switch_cost_all_VM=zeros(1,NumSimulazioni);
channel_cost_all_VM=zeros(1,NumSimulazioni);
tot_cost_all_VM=zeros(1,NumSimulazioni);


f_VM_ultima=zeros(1,NumSimulazioni);


%inizio le simulazioni
P_net_completo=[1:0.25:10000];
%P_net_completo=zeros(1,100);


%Jobs=6 + (10-6).*rand(1,1000); %carico uniformemente distribuito
%Jobs=10.*ones(1,10);  %carico deterministico
%Jobs=ciel(abs(10.*Nrand(0,1)));
%Jobs=random('Poisson',1:NumSimulazioni,1,NumSimulazioni);
%Jobs=ceil(abs(random('Normal',0,4,1,10)));
Jobs=load('workloads8PMR125No2000.mat','L_tot_vect');
Jobs=Jobs.L_tot_vect;
for k=1:NumSimulazioni
    tic
    k
    M=VM(k);
    costo_speed=0;
    costo_switch=0;
    costo_channel=0;
    costo_tot=0;
    
    costo_speed_Nicola=0;
    costo_switch_Nicola=0;
    costo_channel_Nicola=0;
    costo_tot_Nicola=0;
    rate_ultima_macchina=0;
    
    
    k_e=0.005;
    C_max=15; %(Mb/s)
    T_tot=5; %(s)
    
    W = ones(1,M);                  %band of each channel (MHz)
    Zeta = 0.5 * ones(1,M);         %channel parameters (mW)
    
    
    f_max=105.*ones(1,M);  %(Mb/s)
    f_zero=zeros(1,M);  % (Mb/s)
    %f_zero=0.2.*f_max;
    %f_zero=f_max;
    P_net=P_net_completo(1:M); %(mW)
    %P_net=1.*ones(1,M);   %NB.! ORDINARE LE MACCHINE VIRTUALI IN MODO DA AVERE LE POTENZE IN ORDINE CRESCENTE
    %(unicamente perchè semplifica la scritura del
    %software - controindicazione: il codicie funziona bene indipendentemente)
    Th=2.*P_net./C_max;   %Soglia di ibernazione macchina; (tutte le VM per cui mu_opt<=Th vengono ibernate con coefficiente di
    % ibernazione alpha_zero)
    
    E_max=60.*ones(1,M);  % (mJ)
    omega=1.*ones(1,M);
    Delta=0.1.*ones(1,M); % (s)
    Delta_max=max(Delta);
    L_b=zeros(1,M);  % NB. il codice funziona bene per L_b=0, altrimenti aggiungere modifica su f_opt.
    
    temp=2.*k_e+(2.*E_max.*omega./(f_max.^2));
    alpha_zero=2.*k_e./temp;
    alpha_mu=Delta./temp;
    
speed_cost_all_WL=zeros(1,length(Jobs));
switch_cost_all_WL=zeros(1,length(Jobs));
channel_cost_all_WL=zeros(1,length(Jobs));
tot_cost_all_WL=zeros(1,length(Jobs));

    
    
    for num_job=1:length(Jobs)
        L_tot=Jobs(num_job);
        
        
        
        
        %CHECK FEASIBILITY - controllare funzionamento operatori %%%%%%%%%%%%%%%%%%
        
        condizione_feas_carico=(sum(f_max.*Delta-L_b)>=L_tot);
        condizione_feas_tempo=(L_tot<=C_max.*(T_tot-Delta_max)./2);
        condizione_feas_back=min(Delta.*f_max>=L_b);
        
        
        feasibility=condizione_feas_carico && condizione_feas_tempo && condizione_feas_back;
        
        if ~feasibility
            'VM='
            M
            error('Problem unfeasible')
        else
            %'Problem Feasibile!'
        end
        
        
        
        
        
        %SOLUZIONE EQUAZIONE:   sum(L_opt)=L_tot;
        
        %NB. IL PUNTO DI ATTRAVERSAMENTO DELLO ZERO POTREBBE NON ESISTERE. IN QUEL
        %CASO mu_opt E' IL PUNTO DI GRADINO FRA LA ZONA POSITIVA E QUELLA NEGATIVA
        %E delta_mu (CHE E' IL VALORE DELL'EQUAZIONE IN mu) E' DIVERSO DA ZERO.
        
        [mu,delta_mu]= Mu_opt_bisezione(alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
        tol_mu=10^(-2);%tolleranza per decisione macchine in stato limite
        tol_carico_allocato=10^(-2);
        if ((abs(delta_mu)./L_tot)<tol_carico_allocato)%tutto il carico è allocato con un errore relativo inferiore a tol_allocazione
            caso_limite=0;
            %'Caso standard: nessuna VM in stato limite'
        else
            caso_limite=1;
            %'Una o più VM è in stato limite'
            VM_limite=find(abs(Th-mu)<tol_mu);
            if length(VM_limite)==0
                'VM='
                M
                error('nessuna VM in stato limite e carico complessivo L_tot non allocato correttamente')
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%   scheduler ottimo  %%%%%%%%%%%%%%%%%
        
        
        
        
        f_mu=max(mu-2.*P_net./C_max,0);
        f_star=alpha_zero.*f_zero+alpha_mu.*f_mu;
        f_opt=max(0,min(f_star,f_max)); %N.B. --> NEL CASO GENERALE DI L_b>0 IL MAX VA FATTO RISPETTO A L_b/Delta, NON RISPETTO A 0.
        
        rate_ultima_macchina=rate_ultima_macchina+f_opt(M);
        
        
        
        canali_attivi=f_mu>0;
        L_opt=canali_attivi.*(f_opt.*Delta-L_b);
        
        %riallocazione delle eventuali VM sulla soglia (stato limite)
        
        if caso_limite
            L_opt(VM_limite)=0;
            L_allocato=sum(L_opt);
            L_residuo=L_tot-L_allocato;
            for k2=1:length(VM_limite);
                L_opt(VM_limite(k2))=min(L_residuo,f_opt(VM_limite(k2)).*Delta(VM_limite(k2))-L_b(VM_limite(k2)));
                L_residuo=L_residuo-L_opt(VM_limite(k2));
            end
            canali_attivi=(L_opt>0);
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%  costo  %%%%%%%%%%%%%%%%%%
        
        Delta_carico=sum(L_opt)-L_tot;
        if ((abs(Delta_carico)./L_tot)>=tol_carico_allocato)
            'Errore: Carico non correttamente allocato'
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%descrete
        f_precedente=zeros(1,M);
        f_successiva=zeros(1,M);
        x1=zeros(1,M);
        
        Q=6; %numero frequenze discrete;
        f_discrete=zeros(M,Q);
        for conta=1:M
            %f_discrete(conta,:)=[0:f_max(conta)./(Q-1):f_max(conta)];
            f_discrete(conta,:)=[0, 5,50,70,90,105]; % Mhz or Mb/s
        end
        
        f_opt_new=f_opt;
        
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %IDEAL
        costo_speed=costo_speed+sum(((f_opt./f_max).^2).*omega.*E_max);
        costo_switch=costo_switch+sum(k_e.*(f_opt-f_zero).^2);
        costo_channel=costo_channel+sum(2.*P_net.*L_opt./C_max);
        
        %NICOLA
        
        costo_speed_Nicola=costo_speed_Nicola+sum(((((f_precedente/f_max).^2)).*omega.*E_max).*(1-x1)+((((f_successiva/f_max).^2)).*omega.*E_max).*x1);
        costo_switch_Nicola=costo_switch_Nicola+sum((k_e.*(f_precedente-f_zero).^2).*(1-x1)+(k_e.*(f_successiva-f_zero).^2).*x1);
        costo_channel_Nicola=costo_channel_Nicola+((T_tot-Delta(1)).*sum(Zeta .* (2.^(2*L_opt./((T_tot-Delta).*W))-1)));
        %costo_channel_Nicola=costo_channel_Nicola+sum(2.*P_net.*L_opt./C_max);
        
        speed_cost_all_WL(num_job)=sum(((((f_precedente/f_max).^2)).*omega.*E_max).*(1-x1)+((((f_successiva/f_max).^2)).*omega.*E_max).*x1);
        switch_cost_all_WL(num_job)=sum((k_e.*(f_precedente-f_zero).^2).*(1-x1)+(k_e.*(f_successiva-f_zero).^2).*x1);
        channel_cost_all_WL(num_job)=((T_tot-Delta(1)).*sum(Zeta .* (2.^(2*L_opt./((T_tot-Delta).*W))-1)));
        tot_cost_all_WL(num_job)=speed_cost_all_WL(num_job)+switch_cost_all_WL(num_job)+channel_cost_all_WL(num_job);

        
        
        f_zero=f_opt;
        
    end
    
    
    % risultati e deallocazione
    
    
    k
    %IDEAL
    costo_speed_VM(k)=costo_speed./length(Jobs);
    costo_switch_VM(k)=costo_switch./length(Jobs);
    costo_channel_VM(k)=costo_channel./length(Jobs);
    costo_tot_VM(k)=costo_speed_VM(k)+costo_switch_VM(k)+costo_channel_VM(k);
    %NICOLA (COMNET)
    speed_cost_all_VM(k)=costo_speed_Nicola./length(Jobs);
    switch_cost_all_VM(k)=costo_switch_Nicola./length(Jobs);
    channel_cost_all_VM(k)=costo_channel_Nicola./length(Jobs);
    tot_cost_all_VM(k)=speed_cost_all_VM(k)+switch_cost_all_VM(k)+channel_cost_all_VM(k);
    f_VM_ultima(k)=rate_ultima_macchina./length(Jobs);
    
    
    clear f_max f_zero P_net Th  E_max omega Delta L_b temp alpha_zero alpha_mu
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%   istogrammi  %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  figure(k);
    %  bar(f_opt),xlabel('VM');
    %  ylabel('Speed f');
    %  grid on;
    
    toc
end


figure(1)
plot(VM,costo_speed_VM,'--o',VM,speed_cost_all_VM,'--*'),xlabel('VM'),ylabel('Speed Cost');
legend('IDEAL(no DVFS)_type2', 'DVFS_type2');
figure(2)
plot(VM,costo_switch_VM,'--o',VM,switch_cost_all_VM,'--*'),xlabel('VM'),ylabel('Switch Cost');
legend('IDEAL(no DVFS)_type2', 'DVFS_type2');
figure(3)
plot(VM,costo_channel_VM,'--o',VM,channel_cost_all_VM,'--*'),xlabel('VM'),ylabel('Net Cost');
legend('IDEAL(no DVFS)_type2', 'DVFS_type2');
figure(4)
plot(VM,costo_tot_VM,'--o',VM,costo_tot_VM,'--*'),xlabel('VM'),ylabel('Overall Cost');
legend('IDEAL(no DVFS)_type2', 'DVFS_type2');


figure(5)
plot(1:1:length(Jobs),speed_cost_all_WL,'--o'),xlabel('#WL'),ylabel('Speed Cost');
legend('DVFS_type2');
figure(6)
plot(1:1:length(Jobs),switch_cost_all_WL,'--*'),xlabel('#WL'),ylabel('Switch Cost');
legend('DVFS_type2');
figure(7)
plot(1:1:length(Jobs),channel_cost_all_WL,'--*'),xlabel('#WL'),ylabel('Net Cost');
legend('DVFS_type2');
figure(8)
plot(1:1:length(Jobs),tot_cost_all_WL,'--*'),xlabel('#WL'),ylabel('Overall Cost');
legend('DVFS_type2');



figure(9)
plot(VM,f_VM_ultima,'--o'),xlabel('VM'),ylabel('rate last VM');



