
function [ iter_consolidate, L_opt, f_opt, r_opt, m_opt] = Consolidate_operation_baccarelli1(f_opt_ON, f_max, k_e, teta, alpha, g, sigma, E_idle, E_max, T_ON,...
    Delta, f_zero, VM_status_vector, M, L_max, L_tot, N_o, T_tot,W, L_opt_old,mu_old, f_opt_old, r_opt_old, gamma_max, TH, kappa, c, tol_workload_allocated, l_iteratiuons )
mu_iter=zeros(1,l_iteratiuons);
gamma_iter=zeros(1,l_iteratiuons);
V_iter=zeros(1,l_iteratiuons);
L_opt_iter=zeros(M,l_iteratiuons);
r_opt_iter=zeros(M,l_iteratiuons);
f_opt_iter=zeros(M,l_iteratiuons);
L_tild_iter=zeros(M,l_iteratiuons);
L_bar_iter=zeros(M,l_iteratiuons);
psi_hat_iter=zeros(M,l_iteratiuons);
h_hat_iter=zeros(M,l_iteratiuons);
a_hat_iter=zeros(M,l_iteratiuons);
b_hat_iter=zeros(M,l_iteratiuons);

%initilization
mu_iter(1)=mu_old;
f_opt_iter(:,1)=f_opt_old;
L_opt_iter(:,1)=L_opt_old;
L_tild_iter(:,1)=L_opt_old;
r_opt_iter(:,1)=r_opt_old;
psi_hat_iter(:,1)=gamma_max;
gamma_iter(1)=gamma_max;
h_hat_iter(:,1)=gamma_max;
V_iter(1)=0;
a_hat_iter(:,1)=0;
Delta_load=zeros(1,l_iteratiuons);
l=2;
while l<=l_iteratiuons
    y=sum(L_opt_iter(:,l-1))-L_tot;
    
    % calculate mu_iter
    mu_iter(l)=max(0,mu_iter(l-1)-gamma_iter(l-1)*(y));
    
    % calculate f_opt_iter
    [pi_hat_temp, f_differ_result]=pi_hat(f_opt_iter(:,l-1)', f_opt_ON, r_opt_iter(:,l-1)', f_max, k_e, teta, alpha, g, sigma, E_idle, E_max, T_ON,...
        Delta, f_zero, VM_status_vector);
    %f_opt_temp=f_opt_iter(:,l-1)'-h_hat_iter(:,l-1)'.*f_differ_result;
    f_opt_temp=f_opt_iter(:,l-1)'-(kappa/l).*f_differ_result;
    
    f_opt_iter(:,l)=max(0, min(f_max,f_opt_temp));
    
    % calculate L_tild_iter
    L_differ_result=func_Lopt(L_tild_iter(:,l-1)',L_max, N_o, W, T_tot, Delta, g, teta, alpha, mu_iter(l));
    
    %L_tild_temp=L_tild_iter(:,l-1)'-(mu_iter(l)>=TH).*psi_hat_iter(:,l-1)'.*L_differ_result;
    L_tild_temp=L_tild_iter(:,l-1)'-(mu_iter(l)>=TH).*(kappa/l).*L_differ_result;
    L_tild_iter(:,l)=max(0, min(L_max,L_tild_temp));
    
    % calculate L_bar_iter
    
    U_temp=U_1(f_opt_iter(:,l)',f_opt_ON, sigma);
    L_bar_iter(:,l)=Delta.*f_opt_iter(:,l)'-(VM_status_vector'.*(f_opt_iter(:,l)'.*T_ON).*U_temp);
    
    % calculate L_opt_iter
    
    L_opt_iter(:,l)=min(L_bar_iter(:,l),((mu_iter(l)<=TH)*0+(mu_iter(l)>TH).*L_tild_iter(:,l)));
    %L_opt_iter(:,l)=min(L_bar_iter(:,l),((mu_iter(l)>TH).*L_tild_iter(:,l)));
    
    
    % calculate r_opt_iter
    r_differ_result=(-1)*func_Lopt(L_tild_iter(:,l)',L_max, N_o, W, T_tot, Delta, g, teta, alpha, mu_iter(l));
    r_opt_iter(:,l)= max(0,r_differ_result);
    
    % calculate gamma_iter
    
    gamma_iter(l)=max(0,min(gamma_max, gamma_iter(l-1)-c.*(kappa/l)*(y)));
    
  
    
    iter_consolidate=l;
    
    Delta_load(l)=sum(L_opt_iter(:,l))-L_tot;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%hemogenous case
    
   % cond_car_all = (abs(Delta_load(l)/L_tot)<=tol_workload_allocated);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%hetereogenous case

    cond_car_all = (abs(Delta_load(l))<=tol_workload_allocated);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cond_car_all
        m_opt=mu_iter(l);
        L_opt=L_opt_iter(:,l);
        f_opt=f_opt_iter(:,l);
        r_opt=r_opt_iter(:,l);
        break
    else
        %fprintf('incorrect allocation for workload n: %d at iteration %d\n',L_tot,l);
        m_opt=mu_iter(l);
        L_opt=L_opt_iter(:,l);
        f_opt=f_opt_iter(:,l);
        r_opt=r_opt_iter(:,l);
    end
      l=l+1;

end


end



                    