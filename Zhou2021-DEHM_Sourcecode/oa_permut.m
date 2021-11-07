function A = oa_permut(q,n,j)
% Generate an Orthogonal Array using simple permutation method. The
% algorithm was provided in the paper Leung et al.
%
% Calling A=OA_PERMUT(Q,N,J)
%
% Q = the number of levels
% N = the number of factors or columns
% J is selected such that number of rows is equal to Q^J and
% N = (Q^J - 1)/(Q-1)
%
% Tested on oa_permut(3,4,2) given in paper
% Tested on Taguchi L64 array oa_permut(2,63,6)

% Reference: Leung, Y.-W.; Yuping Wang, "An orthogonal genetic algorithm 
% with quantization % for global numerical optimization," Evolutionary 
% Computation, IEEE Transactions on , vol.5, % no.1, pp.41,53, Feb 2001, 
% doi: 10.1109/4235.910464

% Natasha Y Jeppu, natasha.jeppu@gmail.com

if n ~= (q^j-1)/(q-1)
    disp('Does not satisfy criteria ..');
    A=[];
    return
end
row=q^j;
col=(q^j-1)/(q-1);
A=zeros(row,col);

% Compute the basic columns
for k = 1:j
    J=((q^(k-1)-1)/(q-1))+1;
    for i = 1:q^j
        A(i,J)=floor(((i-1)/(q^(j-k)))); % I have to use floor to get the correct result
     end  
end

% Compute the non basic columns
for k = 2:j
    J=((q^(k-1)-1)/(q-1))+1;
    for s = 1:J-1
        for t = 1:q-1
            x=J+(s-1)*(q-1)+t;
            A(:,x)=mod(A(:,s)*t+A(:,J),q);
        end
    end
end
A=mod(A,q);