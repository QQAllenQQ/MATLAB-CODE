function [power_opt,CNR_eff] = Multi_Ns_WFPA(P_total_user,F_HBF,W_HBF,G,U,Ns_alloc,N_sub,N0)
err = 1e-4;
CNR_eff = zeros(max(Ns_alloc),N_sub,U);
power_opt = zeros(max(Ns_alloc),N_sub,U);
for u=1:1:U
    for ns=1:1:Ns_alloc(u)
        for k=1:1:N_sub
            CNR_eff(ns,k,u) = abs(W_HBF(:,ns,u,k)'*G(:,:,k,u)*F_HBF(:,ns,u,k))^2/(N0*norm(W_HBF(:,ns,u,k)')^2);
        end
    end % end of data stream index for u-th user
       %% Tseng
    % Initialization of water filling parameter
    % Set the range of wlevel solution
    ub = max(max(1./CNR_eff(1:Ns_alloc(u),:,u))) + P_total_user;
    lb = 0;
    wlevel = (ub+lb)/2;
    power = sum(sum(max(wlevel-1./CNR_eff(1:Ns_alloc(u),:,u),0)));
    % Bisection method
    while (abs(P_total_user-power) > err) && (wlevel~=ub) && (wlevel~=lb)
        if power >= P_total_user
            ub = wlevel;
        else
            lb = wlevel;
        end
        wlevel = (ub+lb)/2;
        power = sum(sum(max(wlevel-1./CNR_eff(1:Ns_alloc(u),:,u),0)));
    end
    if (wlevel==ub) || (wlevel==lb)
        u
        power_opt(1:Ns_alloc(u),:,u) = P_total_user/N_sub/Ns_alloc(u)*ones(Ns_alloc(u),N_sub);
    else
        power_opt(1:Ns_alloc(u),:,u) = max(wlevel-1./CNR_eff(1:Ns_alloc(u),:,u),0);    
    end 
end % end of user index 
