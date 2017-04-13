function  [Gcl,E,EGS,Psi,Psi_GS,NSz_GS,Problem_mat,nd,ndup,nddown,nc,ncup,ncdown,D,an_m,bn2_m,dplusd,an_p,bn2_p,ddplus] = ED_Green_final(wn,ed,U,ee,VV,Ns,C_ind,table,indice_sector,H_non_zero_ele,spar)       
%Louis-Francois Arsenault Columbia University 2015
%This function calculates all the properties of the ED solution

        %disp('Calculation of the Ground state energy and eigenvector')
        [E,EGS,Psi,Psi_GS,NSz_GS,Problem_mat] = ED_Ns_AIM_final(ed,U,ee,VV,Ns,C_ind,table,indice_sector,H_non_zero_ele,spar);
        [nd,ndup,nddown,nc,ncup,ncdown,D] = Get_nd_D_final(Psi_GS,NSz_GS,C_ind,table,Ns);
    
        %disp('Calculation of the continued fraction expansions')
        [an_m,bn2_m,dplusd] = cont_fract_coeff_G_dGS_final(Psi_GS,NSz_GS,C_ind,table,H_non_zero_ele,ed,U,ee,VV,Ns);
        [an_p,bn2_p,ddplus] = cont_fract_coeff_G_dpluGS_final(Psi_GS,NSz_GS,C_ind,table,H_non_zero_ele,ed,U,ee,VV,Ns);
       
        Gcl = zeros(1,length(wn));
        for r_deg = 1:size(an_m,1)
            Gcl = Gcl + Green_cont_frac_backward_final(an_m(r_deg,:),bn2_m(r_deg,2:size(bn2_m,2)),dplusd(r_deg),EGS(r_deg),-1,i*wn)+Green_cont_frac_backward_final(an_p(r_deg,:),bn2_p(r_deg,2:size(bn2_p,2)),ddplus(r_deg),EGS(r_deg),1,i*wn);
        end
        Gcl = (1/size(an_m,1))*Gcl;
        
        for r_deg1 = 1:size(an_m,1)
            bn2_m(r_deg1,1) = [dplusd(r_deg1)];
            bn2_p(r_deg1,1) = [ddplus(r_deg1)];
        end
        
end
