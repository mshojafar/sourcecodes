function [result, deffer_result] = pi_hat( f_opt, f_opt_ON, r_opt, f_max, k_e, teta, alpha, g, sigma, E_idle, E_max, T_ON,...
                                            Delta, f_zero, VM_status_vector)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%thereshold=10^(-2);
psi_temp=psi_func(f_opt./f_max);

U_temp=U_1(f_opt,f_opt_ON, sigma);
div_tanh_temp=div_tanh( f_opt,f_opt_ON,sigma);
result1=(alpha*teta).*[(E_idle./2.*div_tanh_temp)+((E_max-E_idle)./f_max).*psi_temp+2*k_e.*f_opt];
result2= 2*g.*[max(0,(f_opt-f_max))-max(0,(-1).*f_opt)];
result3=VM_status_vector'.*(r_opt.*T_ON).*U_temp;

result=result1+result2+result3;

deffer_result=result-(alpha*teta*2*k_e.*f_zero+r_opt.*Delta);

end

