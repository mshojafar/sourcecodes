%calcolo moltiplicatori con metodo bisezione

function [mu,delta_mu]= Mu_opt_bisezione(alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot)


%NB.
%LA FUNZIONE DI CUI CALCOLARE LO ZERO è ASSUNTA MONOTONA CRESCENTE
%CON VALORE IN ZERO NEGATIVO E ALL'INFINITO POSITIVO
%SI ASSUME CHE VENGA CHIAMATA SOLO NEI CASI IN CUI LO ZERO ESISTE FINITO.
%ALTRIMENTI VA IN LOOP!

% definisco gli estremi in cui cercare 
a_mu = 0;% per mu=0 si ha sum(L)<L_tot 
%b_mu = 2.*10^3;  %<====  NB. SCEGLIERE CON ATTENZIONE IN BASE AI PARAMETRI DI SISTEMA!
b_mu=max(((1./alpha_mu).*(f_max-alpha_zero.*f_zero))+2.*P_net./C_max); % <== per questo mu si ha f=f_max quindi sum(L)>L_tot
%b_mu = 2.*10^0;


toll = 10^-6;

fa=delta_carico(a_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
fb=delta_carico(b_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);

% verifico che la funzione abbia uno zero  
while (fa.*fb) > 0 
    b_mu=2.*b_mu;
    fb=delta_carico(b_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
    %fprintf('%d\n',t)
    error('estremi iniziali ricerca mu_opt errato')
  
end


iterazioni = ceil(log2(b_mu-a_mu)-log2(toll));  
%     fprintf('%d\n......iteraion',iterazioni);
%     fprintf('%d\n......b_mu',b_mu);
%     fprintf('%d\n......a_mu',a_mu);
%     fprintf('%d\n......toll',toll);

   for i = 1: iterazioni
       c_mu=(a_mu+b_mu)/2;% punto medio
       fc = delta_carico(c_mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
       if abs(fc)<toll
           break
       end
        if abs(b_mu-a_mu)< toll
            break
        end      
%         tolf = toll.*(abs((fb-fa)./(b_mu-a_mu)));     %ELIMINATO - NEL CASO
%                                                       DI GRADINO FUORVIANTE
%         if abs(fc)<=tolf
%             'Attenzione Probabile Gradino'
%             i
%             break
%         end
        if (fa.*fc)<0 %la soluzione è fra a e il punto medio
            b_mu=c_mu;     %tengo la met sinistra dell'intervallo
            fb = fc;
        else         % altrimenti tengo la metà destra
            a_mu=c_mu;
            fa = fc;
        end
   end
    mu=c_mu;
    delta_mu=fc;
