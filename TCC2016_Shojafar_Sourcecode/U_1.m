function [ result ] = U_1( x,a, sigma) %%%%== U_{-1}
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    result=0.5.*(1+tanh((x-a)./sigma));

end

