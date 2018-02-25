%FUNZIONE CHE DEFINISCE L'EQUAZIONE DA RISOLVERE PER IL CALCOLO DEL Mu
%ottimo (equazione di carico su L_tot)




function y=delta_carico(mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot)

f_mu=max(mu-2.*P_net./C_max,0); 
f_star=alpha_zero.*f_zero+alpha_mu.*f_mu;
f_opt=max(0,min(f_star,f_max));

canali_attivi=(f_mu>0); 
L_opt=canali_attivi.*(f_opt.*Delta-L_b);
y=sum(L_opt)-L_tot;  % N.B. PER MU=0 Y<O; PER MU-->INF Y>0  
 