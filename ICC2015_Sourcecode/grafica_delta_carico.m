function y=grafica_delta_carico(mu,alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot)

y=zeros(1,length(mu));
for k=1:length(mu)
    y(k)=delta_carico(mu(k),alpha_zero,alpha_mu,P_net,C_max,f_zero,f_max,Delta,L_b,L_tot);
end