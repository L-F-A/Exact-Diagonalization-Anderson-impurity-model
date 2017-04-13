function [an,bn2,dplusd] = cont_fract_coeff_G_dGS_final(Psi_GS,NSz_GS,C_ind,table,H_non_zero_ele,ed,U,ee,V,Ns)

nbr_coeff_cont = 800;

for r_deg = 1:size(NSz_GS,1)
    Psi_GS_vec = Psi_GS{1,r_deg};
    N_elec_GS = NSz_GS(r_deg,1);
    Sz_GS = NSz_GS(r_deg,2);
    
    
    %We isolate in C the states that are part of the GS
    C_GS_inter = C_ind{1,N_elec_GS+1};
    
    indice_GS = find(C_GS_inter(:,2) == Sz_GS);
    C_GS = C_GS_inter(indice_GS,:);
    
    indice_f0 = find(bitget(C_GS(:,1),1) == 1);
    f01 = Psi_GS_vec(indice_f0,1);
    C_f0 = C_GS(indice_f0,:);
    C_f0_bin = bitset(table(C_f0(:,1)+1,4),2*Ns,0);
    C_f0(:,1) = table(C_f0_bin+1,5);
    C_f0(:,2) = C_f0(:,2) - 1;
    
    C_Lanczos = C_ind{1,N_elec_GS};
    indice_Lanczos_basis = find(C_ind{1,N_elec_GS}(:,2)==(Sz_GS-1));
    C_Lanczos = C_Lanczos(indice_Lanczos_basis,:);
    
    Sz_N = C_ind{1,N_elec_GS}(:,2);
    change_spin = find(diff(Sz_N)~=0)';
    nbr_change = length(change_spin);
    if nbr_change ~=0
      hh = [change_spin(1:length(change_spin)) change_spin(length(change_spin))+1];
      pos_spin_space = find(Sz_N(hh) == Sz_GS-1);
    else
      
    end
    
    
    [lig_f0,col_f0] = size(C_f0);
    [lig_Lan,col_Lan] = size(C_Lanczos);
    

    state_ind_Lan = C_Lanczos(:,1);
    state_ind_f0 = C_f0(:,1);

    f0 = zeros(lig_Lan,1);
    
    if length(f0) == 0
        an(r_deg,1) = 0;
        bn2(r_deg,1) = 0;
        dplusd(r_deg) = 0;
        fprintf('Size of Lanczos basis for sector d_up[%d\t%d]=[%d\t%d]: %d\n\n',N_elec_GS,Sz_GS,N_elec_GS-1,Sz_GS-1,length(indice_Lanczos_basis))
        
    else
    
    for r = 1:length(state_ind_f0) 
        f0(find(state_ind_Lan == state_ind_f0(r))) = f01(r);
    end
    
    fprintf('Size of Lanczos basis for sector d_up[%d\t%d]=[%d\t%d]: %d\n\n',N_elec_GS,Sz_GS,N_elec_GS-1,Sz_GS-1,length(indice_Lanczos_basis))
    
    
                
                nd_states = bitget(table(state_ind_Lan(:,1)+1,4),2*Ns)+bitget(table(state_ind_Lan(:,1)+1,4),2*Ns-1);
                D_states = bitget(table(state_ind_Lan(:,1)+1,4),2*Ns).*bitget(table(state_ind_Lan(:,1)+1,4),2*Ns-1);
                nc_states = zeros(lig_Lan,Ns-1);
                
                nc_states(lig_Lan,Ns-1) = 0;
                for oo = 1:(Ns-1)
                    nc_states(:,oo) = bitget(table(state_ind_Lan(:,1)+1,4),2*Ns-2*oo)+bitget(table(state_ind_Lan(:,1)+1,4),2*Ns-2*oo-1);
                end
                
                diag_H = ed*nd_states + U*D_states + nc_states*ee';
                L_DH = length(diag_H);
                ind_diag = 1:L_DH;
                clear nc_states;
                
                if lig_Lan > 1
                    
                 lig_col_hij = H_non_zero_ele{1,N_elec_GS}{1,pos_spin_space};    
                 lig = lig_col_hij(2,:);
                 col = lig_col_hij(3,:);
                 hij = lig_col_hij(1,:);

                 lig = [lig_col_hij(2,:) ind_diag];
                 col = [lig_col_hij(3,:) ind_diag];
                 hij = sign(lig_col_hij(1,:)).*V(abs(lig_col_hij(1,:)));
                 hij = [hij 0.5*diag_H'];
                 
                 H = sparse(lig,col,hij,L_DH,L_DH);
                 H =  (H+H');

                else
                    H = diag_H;
                end
    
        
        f0_s_f0 = f0'*f0;
        dplusd(r_deg) = f0_s_f0; 
        f0 = f0./sqrt(f0_s_f0);
            
        a0 = (f0'*H*f0);
        b02 = 0;
    
        f1 = H*f0 - a0*f0;
       if f1 == 0
           an(r_deg,1) = a0;
           bn2(r_deg,1) = 0;
           bn2(r_deg,2) = 0;
       else
        b12 = (f1'*f1);
        f1 = f1./sqrt(b12);
        a1 = (f1'*H*f1);
    
    
        an(r_deg,1) = a0;
        an(r_deg,2) = a1;
        bn2(r_deg,1) = 0;
        bn2(r_deg,2) = b12;
        boucle_inf = 1;
        pos = 3;
        
        while boucle_inf == 1
        %pos
            f = H*f1 - a1*f1 - sqrt(b12)*f0;
      
            b = f'*f;
        f = f./sqrt(b);
        a = (f'*H*f);
      
      
        an(r_deg,pos) = a;
        bn2(r_deg,pos) = b;
        pos = pos+1;
        f0 = f1;
        f1 = f;
        a1 = a;
        b12 = b;
        if (b <= 1e-14) || pos == lig_Lan || pos > nbr_coeff_cont
            boucle_inf = 0;
        end
        end
    end
    end
    clear C_f0 nc_states;
end   
    
