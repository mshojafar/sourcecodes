
%VERSIONE DEL TRACKING CON PUNTO INIZIALE DI RICERCA PER MU.

function [conv,iter,mu_opt,ni_opt,f_opt,L_opt,alpha_step_finale,V_step_finale,mu_vect]= tracking3(L_tot,L_b,alpha_zero,alpha_mu,chi,W,T_tot,Delta,f_zero,f_max,M,Th,mu_iniz,ni_iniz,f_iniz,L_iniz,alpha_step_iniz,V_step_iniz)



epsilon_ni=0;  %<== N.B. shift numerico introdotto su ni;

toll = 10^(-3);%10^(-2);
num_iterazioni=100000;
conv=0; %viene settato a 1 se il tracking arriva a convergenza.

ni=zeros(num_iterazioni,M);
ni(1,:)=ni_iniz;
f=zeros(num_iterazioni,M);
f(1,:)=f_iniz;
L=zeros(num_iterazioni,M);
L(1,:)=L_iniz;
mu=zeros(num_iterazioni,1);
mu(1)=mu_iniz; %<==  importante, inizializzazione del parametro mu.


beta=0.04;%0.01;%0.7;%0.1;%10^(-3);%10^(-4);
gamma=0.6;%0.4;
alpha_step=beta.*ones(num_iterazioni,1);
%alpha_step(1)=alpha_step_iniz;
V_step=zeros(num_iterazioni,1);
%V_step(1)=V_step_iniz;

gap=zeros(num_iterazioni,1);


for k=2:num_iterazioni
    
    gap(k-1)=sum(L(k-1,:))-L_tot;
    if abs(gap(k-1)./L_tot)<toll   %errore relativo sotto la tolleranza
        conv=1;
        break
    end
    mu(k)=mu(k-1)-alpha_step(k-1).*gap(k-1);
    mu(k)=max(0,mu(k));
    
    ni_star=mu(k)-2.*log(2).*(chi./W).*2.^(2.*L(k-1,:)./(W.*(T_tot-Delta)));
    ni(k,:)=max(0,ni_star);
    
    f_star=alpha_zero.*f_zero+alpha_mu.*ni(k,:);
    f(k,:)=max(0,min(f_star,f_max)); %N.B. --> NEL CASO GENERALE DI L_b>0 IL MAX VA FATTO RISPETTO A L_b/Delta, NON RISPETTO A 0.
    
    canali_attivi=ni(k,:)>epsilon_ni;  %<== N.B. shift numerico introdotto su ni;;
    %L(k,:)=canali_attivi.*(f(k,:).*Delta-L_b)+(~canali_attivi).*max(0,(((T_tot-Delta).*W./2).*log(mu(k)./Th)./log(2)));
    %L(k,:)=max(0,(((T_tot-Delta).*W./2).*log(mu(k)./Th)./log(2)));
    L(k,:)=min((f(k,:).*Delta-L_b),max(0,(((T_tot-Delta).*W./2).*log(mu(k)./Th)./log(2))));
    
    
    
    
    V_step(k)=(1-alpha_step(k-1)).*V_step(k-1)-gap(k-1);
    alpha_step(k)=max(0,(min(beta,alpha_step(k-1)-gamma.*V_step(k-1).*gap(k-1))));
    %alpha_step(k)=(min(beta,alpha_step(k-1)-gamma.*V_step(k-1).*gap(k-1)));
    
end



%  V_step
%  alpha_step
% gap
% mu
% ni
% f
% L


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GRAFICI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%
% 
% figure(num_iterazioni+1)
% plot(gap,'--o'),xlabel('iterazioni'),ylabel('gap');
% 
% figure(num_iterazioni+2)
% plot(mu,'--o'),xlabel('iterazioni'),ylabel('mu');
% 
% figure(num_iterazioni+3)
% plot(ni,'--o'),xlabel('iterazioni'),ylabel('ni');
% 
% figure(num_iterazioni+4)
% plot(f,'--o'),xlabel('iterazioni'),ylabel('f');
% 
% figure(num_iterazioni+5)
% plot(L,'--o'),xlabel('iterazioni'),ylabel('L');
% 
% figure(num_iterazioni+6)
% plot(V_step,'--o'),xlabel('iterazioni'),ylabel('V_step');
% 
% figure(num_iterazioni+7)
% plot(alpha_step,'--o'),xlabel('iterazioni'),ylabel('alpha_step');




%  I VALORI SONO CORRETTI SOLTANTO SE CONV=1
iter=k-1; %numero iterazioni per convergenza
mu_opt=mu(k-1);
mu_vect=mu(1:k-1);% < == N.B. Restituisce in uscita il vettoredi tracking
ni_opt=ni(k-1,:);
f_opt=f(k-1,:);
L_opt=L(k-1,:);

alpha_step_finale=alpha_step(k-1);
V_step_finale=V_step(k-1);
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% 
% 
% 
% %NB.
% %LA FUNZIONE DI CUI CALCOLARE LO ZERO è ASSUNTA MONOTONA CRESCENTE
% %CON VALORE IN ZERO NEGATIVO E ALL'INFINITO POSITIVO
% %SI ASSUME CHE VENGA CHIAMATA SOLO NEI CASI IN CUI LO ZERO ESISTE FINITO.
% %ALTRIMENTI VA IN LOOP!
% 
% % definisco gli estremi in cui cercare 
% a_mu = 0;% per mu=0 si ha sum(L)<L_tot 
% %b_mu = 2.*10^3;  %<====  NB. SCEGLIERE CON ATTENZIONE IN BASE AI PARAMETRI DI SISTEMA!
% b_mu=max(((1./alpha_mu).*(f_max-alpha_zero.*f_zero))+2.*P_net./C_max); % <== per questo mu si ha f=f_max quindi sum(L)>L_tot
% %b_mu = 2.*10^0;
% 
% 
% toll = 10^-6;
% 
% fa=delta_carico(a_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
% fb=delta_carico(b_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
% 
% % verifico che la funzione abbia uno zero  
% while (fa.*fb) > 0 
%     b_mu=2.*b_mu;
%     fb=delta_carico(b_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
%     %fprintf('%d\n',t)
%     error('estremi iniziali ricerca mu_opt errato')
%   
% end
% 
% 
% iterazioni = ceil(log2(b_mu-a_mu)-log2(toll));  
%     
% 
%    for i = 1: iterazioni
%        c_mu=(a_mu+b_mu)/2;% punto medio
%        fc = delta_carico(c_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
%        if abs(fc)<toll
%            break
%        end
%         if abs(b_mu-a_mu)< toll
%             break
%         end      
% %         tolf = toll.*(abs((fb-fa)./(b_mu-a_mu)));     %ELIMINATO - NEL CASO
% %                                                       DI GRADINO FUORVIANTE
% %         if abs(fc)<=tolf
% %             'Attenzione Probabile Gradino'
% %             i
% %             break
% %         end
%         if (fa.*fc)<0 %la soluzione è fra a e il punto medio
%             b_mu=c_mu;     %tengo la met sinistra dell'intervallo
%             fb = fc;
%         else         % altrimenti tengo la metà destra
%             a_mu=c_mu;
%             fa = fc;
%         end
%    end
%     mu=c_mu;
%     delta_mu=fc;
