function [E,EGS,Psi,Psi_GS,NSz_GS,Problem_mat] = ED_Ns_AIM_final(ed,U,ee,V,Ns,C_ind,table,indice_sector,H_non_zero_ele,spar)
     
%        opt.DISP = 0; 
%        warning off;

%opts.p = 2;
    file_name_problem = 'No convergence of eigenvalues after 5 attempts for sector [N Sz] = [%d %d]';
    Problem_mat = [0 0 0];
    prob = 1;
    if spar == 0
        nb_eigen = 1;
        OPTS.tol = 1e-30;%eps;
    end

    f = size(C_ind);
    Nb_sector = f(2);
    %fprintf('Number of sectors in particle number Nb_sector = %d\n\n',Nb_sector)
     
    EGS = 0;
    E = cell(f);
    Psi = cell(f);
     
    temps = tic;  
     for r = 1:Nb_sector
        
        
        Sz_N = C_ind{1,r}(:,2);
        change_spin = find(diff(Sz_N)~=0)';
        indice_change_spin = change_spin+1;
        split_spin_sec = [0 change_spin];
        nbr_change = length(indice_change_spin);
        nbr_spin_diff = nbr_change+1;
        
        if nbr_change ~=0
            hh = [change_spin(1:length(change_spin)) change_spin(length(change_spin))+1];
        end
        
        if nbr_change == 0
         tstart1 = tic;   
            %fprintf('\n\n Calculating sector N = %d and Sz = %d\n',r-1,Sz_N)
            %fprintf('This sector has 1 state\n')
             bin_num = table(C_ind{1,r}(1)+1,4);
             CC = cast(dec2bin(bin_num,2*Ns),'double')-48;
             CC(1,2*Ns+1) = Sz_N;
             E{1,r} = Mat_element_AIM_final(CC,CC,Ns,ed,U,ee,V);
             Psi{1,r} = 1;
             %fprintf('The ground state of the sector is %f\n',E{1,r})
             if E{1,r} < EGS
                 EGS = E{1,r};
                 Psi_GS{1,1} = Psi{1,r};
                 %Psi_GS = Psi{1,r};
                 NSz_GS = [r-1 Sz_N];
             end
        %toc(tstart1)    
           
        else
            E{1,r} = cell(1,nbr_change);
            Psi{1,r} = cell(1,nbr_change);
            
            for h = 1:nbr_change+1;
            tstart1 = tic ;
                
                length_spins = length(find(Sz_N == Sz_N(hh(h))));
                %fprintf('\n\nCalculating sector N = %d and Sz = %d\n',r-1,Sz_N(hh(h)))
                %fprintf('This sector has %d states\n',length_spins)
                
                
                states_ket_ind = C_ind{1,r}(1+split_spin_sec(h):length_spins+split_spin_sec(h),:);
                [lig_states,col_states] = size(states_ket_ind);
                
                nd_states = bitget(table(states_ket_ind(:,1)+1,4),2*Ns)+bitget(table(states_ket_ind(:,1)+1,4),2*Ns-1);
                
                D_states = bitget(table(states_ket_ind(:,1)+1,4),2*Ns).*bitget(table(states_ket_ind(:,1)+1,4),2*Ns-1);
                nc_states = zeros(lig_states,Ns-1);
                
                nc_states(lig_states,Ns-1) = 0;
                for oo = 1:(Ns-1)
                    nc_states(:,oo) = bitget(table(states_ket_ind(:,1)+1,4),2*Ns-2*oo)+bitget(table(states_ket_ind(:,1)+1,4),2*Ns-2*oo-1);
                end
                
                diag_H = ed*nd_states + U*D_states + nc_states*ee';
                L_DH = length(diag_H);
                ind_diag = 1:L_DH;
                clear nc_states; 
                if length_spins > 1
                    
                 lig_col_hij = H_non_zero_ele{1,r}{1,h};    
                 lig = lig_col_hij(2,:);
                 col = lig_col_hij(3,:);
                 hij = lig_col_hij(1,:);

                 lig = [lig_col_hij(2,:) ind_diag];
                 col = [lig_col_hij(3,:) ind_diag];
                 hij = sign(lig_col_hij(1,:)).*V(abs(lig_col_hij(1,:)));
                 hij = [hij 0.5*diag_H'];
                 
                 H = sparse(lig,col,hij,L_DH,L_DH);
                 H =  (H+H');
                 %full(H)
                 

                   % if (spar == 0) && (length_spins > 1000)
                        %disp('yup')
                        %[eigen_val,eigen_vec] = eigifp(H,opt);
                        r_convergence = 0;
                        r_count = 1;
                        while r_convergence == 0;
                            [eigen_vec,eigen_val,Flag] = eigs(H,[],nb_eigen,'SA',OPTS);
                            if Flag ~= 0 && r_count ~= 5
                                Problem_mat(prob,1) = r-1;
                                Problem_mat(prob,2) = Sz_N(hh(h));
                                Problem_mat(prob,3) = Flag;
                                prob = prob + 1;
                            elseif Flag ~= 0 && r_count == 5
                                name_problem = sprintf(file_name_problem,r-1,Sz_N(hh(h)));
                                disp(name_problem)
                                Problem_mat(prob,1) = r-1;
                                Problem_mat(prob,2) = Sz_N(hh(h));
                                Problem_mat(prob,3) = Flag;
                                prob = prob + 1;
                                r_convergence = 1;
                            else
                                r_convergence = 1;
                            end
                            r_count = r_count + 1;
                        end
                    %else
                     %   [eigen_vec,eigen_val] = eig(full(H));
                    %end
                    
                    E{1,r}{1,h} = diag(eigen_val);
                    Psi{1,r}{1,h} = eigen_vec;
                    EGS_sec = min(E{1,r}{1,h}); 
                    %fprintf('The ground state of the sector is %f\n',EGS_sec)
                    if EGS_sec < EGS
                        [EGS,L] = min(E{1,r}{1,h});
                        Psi_GS{1,1} = Psi{1,r}{1,h}(:,L);
                        %Psi_GS = Psi{1,r}{1,h}(:,L);
                        NSz_GS = [r-1 Sz_N(hh(h))];
                    end
                else
                    E{1,r}{1,h} = diag_H;
                    Psi{1,r}{1,h} = 1;
                    %fprintf('The ground state of the sector is %f\n',E{1,r}{1,h})
                    if E{1,r}{1,h} < EGS
                        EGS = E{1,r}{1,h};
                        Psi_GS{1,1} = Psi{1,r}{1,h};
                        %Psi_GS = Psi{1,r}{1,h};
                        NSz_GS = [r-1 Sz_N(hh(h))];
                    end
                end
             %toc(tstart1)   
            end
        end
     end
     NSz_GSS = NSz_GS;
     %EGS
     %fprintf('\n\n The full diagonalization took:\n');
     %toc(temps)
     r_psi = 2;
     disp('Starting to search for degeneracies');
     if NSz_GS(2) ~= 0
         new_SZ = -NSz_GS(2);
         new_N_deg = NSz_GS(1);
         NSz_GS = [NSz_GS;new_N_deg new_SZ]
         Sz_N_max = C_ind{1,new_N_deg+1}(1,2);
         dddd = length(Sz_N_max:-2:new_SZ);
         Psi_GS{1,r_psi} = Psi{1, new_N_deg+1}{1,dddd};
         EGS(r_psi) = E{1,new_N_deg+1}{1,dddd};
         r_psi = r_psi+1;
     end
     %if (NSz_GSS(1) ~= Ns) && (NSz_GSS(2) ~= 0)
     for dd = 1:2*Ns+1
         next_loop = length(E{1,dd});
         for ddd = 1:next_loop
%              if dd == 1 || dd == 2*Ns+1
%                  diff_deg = abs(E{1,dd}-EGS);
%              else
%                 diff_deg = abs(E{1,dd}{1,ddd}-EGS);
%              end
               N_electron = dd-1;
               Sz_N_max = C_ind{1,dd}(1,2);
               Sz_now_vec = Sz_N_max:-2:-Sz_N_max;
               Sz_now = Sz_now_vec(ddd);
               
               if ( (N_electron == NSz_GSS(1) && Sz_now ~= -NSz_GSS(2)) || (N_electron ~= NSz_GSS(1)) )
                   if dd == 1 || dd == 2*Ns+1
                     %diff_deg = abs(E{1,dd}-EGS);
                     a = E{1,dd};
                   else
                    %diff_deg = abs(E{1,dd}{1,ddd}-EGS);
                    a = E{1,dd}{1,ddd};
                   end
                    b = EGS(1);
                    aa = abs(a-floor(a));
                    bb = abs(b-floor(b));
                    %The tow energies are equal up to how many digits
                    diff_EGS = floor(abs(log10(abs(b-a))));
                    if  diff_EGS >= 12 %(( (floor(1e12*aa)- floor(1e12*bb)) == 0 ) && ( (floor(1e13*aa)- floor(1e13*bb)) ~= 0)) || (( (floor(1e13*aa)- floor(1e13*bb)) == 0 ) && ( (floor(1e14*aa)- floor(1e14*bb)) ~= 0)) || (( (floor(1e14*aa)- floor(1e14*bb)) == 0 ) && ( (floor(1e15*aa)- floor(1e15*bb)) ~= 0))  || ( (floor(1e15*aa)- floor(1e15*bb)) == 0) 
                     %(( (floor(1e11*aa)- floor(1e11*bb)) == 0 ) && ( (floor(1e12*aa)- floor(1e12*bb)) ~= 0)) ||
                         N_deg = dd-1;
                         Sz_N = C_ind{1,dd}(:,2);
                         change_spin = find(diff(Sz_N)~=0)';
                         indice_change_spin = change_spin+1;
                         split_spin_sec = [0 change_spin];
                         nbr_change = length(indice_change_spin);
                         nbr_spin_diff = nbr_change+1;

                         if nbr_change ~=0
                            hh = [change_spin(1:length(change_spin)) change_spin(length(change_spin))+1];
                         end
                         if (N_deg == NSz_GSS(1) && Sz_N(hh(ddd)) ~= NSz_GSS(2) ) || N_deg ~= NSz_GSS(1)
                            NSz_GS = [NSz_GS;N_deg Sz_N(hh(ddd))];
                            Psi_GS{1,r_psi} = Psi{1,dd}{1,ddd};
                            %Psi_GS(1,r_psi) = Psi{1,dd}{1,ddd};
                            EGS(r_psi) = a;
                            r_psi = r_psi+1;
                         end
                    end
              end
         end
     end
     NSz_GS
end
%end
     
%end

