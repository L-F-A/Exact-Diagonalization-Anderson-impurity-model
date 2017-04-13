function  [Gcl,E,EGS,Psi,Psi_GS,NSz_GS,Problem_mat,nd,ndup,nddown,nc,ncup,ncdown,D,an_mqio,bn2_mqio,dplusd,an_pqio,bn2_pqio,ddplus] = ED_Green_final(wn,ed,U,ee,VV,Ns,C_ind,table,indice_sector,H_non_zero_ele,spar)       
%Louis-Francois Arsenault Columbia University 2015
%This function calculates all the properties of the ED solution

        %disp('Calculation of the Ground state energy and eigenvector')
        [E,EGS,Psi,Psi_GS,NSz_GS,Problem_mat] = ED_Ns_AIM_final(ed,U,ee,VV,Ns,C_ind,table,indice_sector,H_non_zero_ele,spar);
        [nd,ndup,nddown,nc,ncup,ncdown,D] = Get_nd_D_final(Psi_GS,NSz_GS,C_ind,table,Ns);
    
        %disp('Calculation of the continued fraction expansions')
        [an_mqio,bn2_mqio,dplusd] = cont_fract_coeff_G_dGS_final(Psi_GS,NSz_GS,C_ind,table,H_non_zero_ele,ed,U,ee,VV,Ns);
        [an_pqio,bn2_pqio,ddplus] = cont_fract_coeff_G_dpluGS_final(Psi_GS,NSz_GS,C_ind,table,H_non_zero_ele,ed,U,ee,VV,Ns);
       
        Gcl = zeros(1,length(wn));
        for r_deg = 1:size(an_mqio,1)
            Gcl = Gcl + Green_cont_frac_backward_final(an_mqio(r_deg,:),bn2_mqio(r_deg,2:size(bn2_mqio,2)),dplusd(r_deg),EGS(r_deg),-1,i*wn)+Green_cont_frac_backward_final(an_pqio(r_deg,:),bn2_pqio(r_deg,2:size(bn2_pqio,2)),ddplus(r_deg),EGS(r_deg),1,i*wn);
        end
        Gcl = (1/size(an_mqio,1))*Gcl;
        
        for r_deg1 = 1:size(an_mqio,1)
            bn2_mqio(r_deg1,1) = [dplusd(r_deg1)];
            bn2_pqio(r_deg1,1) = [ddplus(r_deg1)];
        end
        
end
