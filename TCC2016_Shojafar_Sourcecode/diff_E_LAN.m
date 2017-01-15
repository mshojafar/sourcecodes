function [ result] = diff_E_LAN(L_opt, N_o, W, T_s, Delta)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

result= 2.*(log(2)).*N_o.*exp((2*log(2)*L_opt)./(W*(T_s-Delta)));

%TH=2.*(log(2)).*N_o;


end

