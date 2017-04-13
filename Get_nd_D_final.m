function [nd,ndup,nddown,nc,ncup,ncdown,D] = Get_nd_D_final(Psi_GS,NSz_GS,C_ind,table,Ns)

  for r_deg = 1:size(NSz_GS,1);
    %r_deg
    Psi_GS_vec = Psi_GS{1,r_deg};  
    N_elec_GS = NSz_GS(r_deg,1);
    Sz_GS = NSz_GS(r_deg,2);
    
    
    %We isolate in C the states that are part of the GS
    %C_GS_inter = C{1,N_elec_GS+1};
    C_GS_inter = C_ind{1,N_elec_GS+1};
    
    %indice_GS = find(C_GS_inter(:,2*Ns+1) == Sz_GS);
    indice_GS = find(C_GS_inter(:,2) == Sz_GS);
    C_GS = C_GS_inter(indice_GS,:);
    
%     Sz_N = C_ind{1,N_elec_GS+1}(:,2);
%     change_spin = find(diff(Sz_N)~=0)';
%     nbr_change = length(change_spin);
%     if nbr_change ~=0
%       hh = [change_spin(1:length(change_spin)) change_spin(length(change_spin))+1];
%       pos_spin_space = find(Sz_N(hh) == Sz_GS+1);
%     end
    
    [lig_GS,col_GS] = size(C_GS);
    state_ind_GS = C_GS(:,1);
    
    nd_states_up = bitget(table(state_ind_GS(:,1)+1,4),2*Ns);
    nd_states_down = bitget(table(state_ind_GS(:,1)+1,4),2*Ns-1);
    nd_states = nd_states_up + nd_states_down;
    %nd_states = bitget(table(state_ind_GS(:,1)+1,4),2*Ns)+bitget(table(state_ind_GS(:,1)+1,4),2*Ns-1);
    D_states = bitget(table(state_ind_GS(:,1)+1,4),2*Ns).*bitget(table(state_ind_GS(:,1)+1,4),2*Ns-1);
%     nd(r_deg) = Psi_GS(:,r_deg)'*(nd_states.*Psi_GS(:,r_deg));
%     ndup(r_deg) = Psi_GS(:,r_deg)'*(nd_states_up.*Psi_GS(:,r_deg));
%     nddown(r_deg) = Psi_GS(:,r_deg)'*(nd_states_down.*Psi_GS(:,r_deg));
%     D(r_deg) = Psi_GS(:,r_deg)'*(D_states.*Psi_GS(:,r_deg));
    nd(r_deg) = Psi_GS_vec'*(nd_states.*Psi_GS_vec);
    ndup(r_deg) = Psi_GS_vec'*(nd_states_up.*Psi_GS_vec);
    nddown(r_deg) = Psi_GS_vec'*(nd_states_down.*Psi_GS_vec);
    D(r_deg) = Psi_GS_vec'*(D_states.*Psi_GS_vec);
    
    
    nc_states_up(lig_GS,Ns-1) = 0;
    nc_states_down(lig_GS,Ns-1) = 0;
    
    for oo = 1:(Ns-1)
       nc_states_up(:,oo) = bitget(table(state_ind_GS(:,1)+1,4),2*Ns-2*oo);
       nc_states_down(:,oo) = bitget(table(state_ind_GS(:,1)+1,4),2*Ns-2*oo-1);
       nc_states(:,oo) = nc_states_up(:,oo)+nc_states_down(:,oo);
       %nc_states(:,oo) = bitget(table(state_ind_GS(:,1)+1,4),2*Ns-2*oo)+bitget(table(state_ind_GS(:,1)+1,4),2*Ns-2*oo-1);
%        nc_inter(oo) = Psi_GS(:,r_deg)'*(nc_states(:,oo).*Psi_GS(:,r_deg));
%        nc_inter_up(oo) = Psi_GS(:,r_deg)'*(nc_states_up(:,oo).*Psi_GS(:,r_deg));
%        nc_inter_down(oo) = Psi_GS(:,r_deg)'*(nc_states_down(:,oo).*Psi_GS(:,r_deg));
       nc_inter(oo) = Psi_GS_vec'*(nc_states(:,oo).*Psi_GS_vec);
       nc_inter_up(oo) = Psi_GS_vec'*(nc_states_up(:,oo).*Psi_GS_vec);
       nc_inter_down(oo) = Psi_GS_vec'*(nc_states_down(:,oo).*Psi_GS_vec);
    end
    
    nc(r_deg) = sum(nc_inter);
    ncup(r_deg) = sum(nc_inter_up);
    ncdown(r_deg) = sum(nc_inter_down);
    clear nc_states_up nc_states_down nc_states
  end
    