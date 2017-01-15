%calcolo moltiplicatori con metodo bisezione

function [iter,mu,delta_mu,L_opt]= Mu_opt_bisezione1(f_max,Delta,Th,L_tot,W,T_tot)


%NB.
%LA FUNZIONE DI CUI CALCOLARE LO ZERO è ASSUNTA MONOTONA CRESCENTE
%CON VALORE IN ZERO NEGATIVO E ALL'INFINITO POSITIVO
%SI ASSUME CHE VENGA CHIAMATA SOLO NEI CASI IN CUI LO ZERO ESISTE FINITO.
%ALTRIMENTI VA IN LOOP!

% definisco gli estremi in cui cercare 
a_mu = 0;% per mu=0 si ha sum(L)<L_tot 
%b_mu = 2.*10^3;  %<====  NB. SCEGLIERE CON ATTENZIONE IN BASE AI PARAMETRI DI SISTEMA!
%b_mu=max(((1./alpha_mu).*(f_max-alpha_zero.*f_zero))+2.*P_net./C_max); % <== per questo mu si ha f=f_max quindi sum(L)>L_tot
temp1=((T_tot-Delta).*W./2);
temp2=(f_max.*Delta)./temp1;
b_mu=max((2.^temp2).*Th); % <== per questo mu si ha f=f_max quindi sum(L)>L_tot
 
%b_mu = 2.*10^0;


toll = 10^-6;

fa=delta_carico(a_mu,f_max,Delta,Th,L_tot,W,T_tot);
fb=delta_carico(b_mu,f_max,Delta,Th,L_tot,W,T_tot);

% verifico che la funzione abbia uno zero  
while (fa.*fb) > 0 
    b_mu=2.*b_mu;
    fb=delta_carico(b_mu,f_max,Delta,Th,L_tot,W,T_tot);
    %fprintf('%d\n',t)
    error('estremi iniziali ricerca mu_opt errato')
  
end


iterazioni = ceil(log2(b_mu-a_mu)-log2(toll));  
    
iter=iterazioni;
   for i = 1: iterazioni
       c_mu=(a_mu+b_mu)/2;% punto medio
       fc = delta_carico(c_mu,f_max,Delta,Th,L_tot,W,T_tot);
       if abs(fc)<toll
           iter=i;
           break
       end
        if abs(b_mu-a_mu)< toll
            iter=i;
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
    L_opt=min(max(0,(((T_tot-Delta).*W./2).*log(mu./Th)./log(2))), f_max.*Delta);
