function [ result ] = div_tanh( x,a,sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
 
result= 1./(sigma.*(cosh((x-a)./sigma)).^2);

end

