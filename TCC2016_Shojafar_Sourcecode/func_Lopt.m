function [result] = func_Lopt(L_opt, L_max, N_o, W, T_s, Delta, g, teta, alpha, mu_opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


result1= diff_E_LAN(L_opt, N_o, W, T_s, Delta);

result2= 2*g*[max(0,(L_opt-L_max))-max(0,(-1).*L_opt)];

result=teta*alpha*result1+result2-mu_opt;

end

