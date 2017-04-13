function [C_ind,table,indice_sector,H_non_zero_ele] = ED_Ns_generate_final(Ns)
     

    fprintf('\n Generating the basis for Ns = %d\n\n',Ns)
    [C_ind,table,indice_sector] = generate_basis_emb(Ns);
    
    f = size(C_ind);
    Nb_sector = f(2);
    fprintf('Number of sectors in particle number Nb_sector = %d\n\n',Nb_sector)
     
    H_non_zero_ele = cell(f);
     
    temps = tic;  
    for r = 1:Nb_sector
        
        Sz_N = C_ind{1,r}(:,2);
        change_spin = find(diff(Sz_N)~=0)';
        indice_change_spin = change_spin+1;
        split_spin_sec = [0 change_spin];
        nbr_change = length(indice_change_spin);
        %nbr_spin_diff = nbr_change+1;
        
        if nbr_change ~=0
            hh = [change_spin(1:length(change_spin)) change_spin(length(change_spin))+1];
        end
        
        if nbr_change == 0
         tstart1 = tic;   
          %  fprintf('\n\n Calculating sector N = %d and Sz = %d\n',r-1,Sz_N)
           % fprintf('This sector has 1 state\n')
               H_non_zero_ele{1,r} = [0 0 0]; 


       %toc(tstart1)   ; 
           
        else
            
            H_non_zero_ele{1,r} = cell(1,nbr_change);
            
            for h = 1:nbr_change+1;
            tstart1 = tic ;
                
                length_spins = length(find(Sz_N == Sz_N(hh(h))));
            %    fprintf('\n\nCalculating sector N = %d and Sz = %d\n',r-1,Sz_N(hh(h)))
             %   fprintf('This sector has %d states\n',length_spins)
                
                states_ket_ind = C_ind{1,r}(1+split_spin_sec(h):length_spins+split_spin_sec(h),:);
                [lig_states,col_states] = size(states_ket_ind);
                
                
                if length_spins > 1
                    [lig_col_hij] = Mat_element_H_emb(states_ket_ind,table,indice_sector(r-1)+1+split_spin_sec(h),lig_states,Ns);
                    H_non_zero_ele{1,r}{1,h} = lig_col_hij;
                else
                    H_non_zero_ele{1,r}{1,h} = [0 0 0];
                end
             %toc(tstart1) ;  
            end
            
        end
        
       
     end
     %fprintf('\n\n The full diagonalization took:\n');
     %toc(temps);
end

function [lig_col_hij] = Mat_element_H_emb(states_ket_ind,table,start,lig_states,Ns)


  g = 1;

%       lig_col_hij(2,num_non_H) = 0; %= zeros(1,num_non_H);
%       lig_col_hij(3,num_non_H) = 0; %= lig;%zeros(1,num_non_H);
%       lig_col_hij(1,num_non_H) = 0; %= lig;%zeros(1,num_non_H);
 
for ll = 1:lig_states
    
    ket = states_ket_ind(ll,1);
    
    %If the indice is odd, there is one up d electron since only the first
    %position has 1 in II_up
    ndu = bitget(ket,1);
    ndd = bitget(table(ket+1,3),1);
    
  if ( ndu ~= 0 ) || ( ndd ~=0 ) 
    ket_bin = table(ket+1,4);
    for r=2:Ns
            
            
        if ( ndu + 2*bitget(ket_bin,2*Ns-2*r+2) ) == 1
            
            p = sum(bitget(ket_bin,(2*Ns-2*r+3):2*Ns-1));
            phase = -mod(p,2)+mod(p+1,2);
            ket_bin_up = ket_bin;
            ket_bin_up = bitset(ket_bin_up,2*Ns,0);
            ket_bin_up = bitset(ket_bin_up,2*Ns-2*r+2,1);
            new_ind_up = table(ket_bin_up+1,5);
            lig_col_hij(1,g) = phase*(r-1);
            lig_col_hij(2,g) = table(new_ind_up+1,2)-start+1;
            lig_col_hij(3,g) = ll;
            g = g+1;
        end
        if ( ndd + 2*bitget(ket_bin,2*Ns - 2*r + 1) ) == 1                                      
            p = sum(bitget(ket_bin,(2*Ns-2*r+2):2*Ns-2));
            phase = -mod(p,2)+mod(p+1,2);
            ket_bin_down = ket_bin;
            ket_bin_down = bitset(ket_bin_down,2*Ns-1,0);
            ket_bin_down = bitset(ket_bin_down,2*Ns-2*r+1,1);
            new_ind_down = table(ket_bin_down+1,5);
            lig_col_hij(1,g) = phase*(r-1);
            lig_col_hij(2,g) = table(new_ind_down+1,2)-start+1;
            lig_col_hij(3,g) = ll;
            g = g+1;
            
        end
    end
  end
end

end

function [C_ind,table,indice_sector] = generate_basis_emb(Ns)
tic
    x1 = 0:2^(2*Ns)-1;
    x = dec2bin(x1);
    M = cast(x,'double');
    M = M-48;
    [rr,hh] = size(x);
     II_up(1,2*Ns) = 0;
     II_down(1,2*Ns) = 0;
    for rrr = 1:2*Ns
        II_up(rrr) = mod(rrr,2)*2^((rrr+1)/2-1)+mod(rrr+1,2)*2^(Ns+rrr/2-1);
        II_down(rrr) = mod(rrr,2)*2^(Ns+ (rrr-1)/2)+mod(rrr+1,2)*2^(rrr/2-1);
    end
        
    S_vec = (-1).^(0:hh-1);
    
    N = sum(M,2);
    %Sz = M*S_vec';
    
    [N,indice] = sort(N);
    %Sz = Sz(indice);
    M = M(indice,:);
    
    M_II_up = M*II_up';
    x11 = x1(indice)';
    
    [M_II_up_sort,ind_s11] = sort(M_II_up);
    [x11_sort,ind_s1] = sort(x11);
    
     
    sector = 0:2*Ns;
    Nb_sectors = length(sector);
     
    C_ind = cell(1,Nb_sectors);
    
    p = find(diff(N) == 1);
    debut = p+1;
    
     
    C_ind{1,1}(1,1) = 0;
    C_ind{1,1}(1,2) = 0;
    
    C_ind{1,Nb_sectors}(1,1) = sum(II_up);
    C_ind{1,Nb_sectors}(1,2) = 0;
    M1(size(M,1),2*Ns) = 0;
    M1(1,:) = M(1,:);
    M1(4^Ns,:) = M(rr,:);
    Linter = 1;
    indice_sector(1) = Linter;
    g = 2;
    
    
    for r = 2:Nb_sectors-1
         
         C = M(debut(r-1):p(r),:);
         C_ind{1,r} = M_II_up(debut(r-1):p(r));
         C_ind{1,r}(:,2) = M(debut(r-1):p(r),:)*S_vec';
         [ordre,ind] = sort(C_ind{1,r}(:,2),'descend');
         C = C(ind,:);
         C_ind{1,r} = C_ind{1,r}(ind,:);
         
         D = size(M(debut(r-1):p(r),:));
         M1(Linter+1:Linter+D(1),:) = C;
         Linter = Linter+D(1);
         indice_sector(g) = Linter;
         g = g+1;
    end
    
    table_inter(4^Ns,5) = 0;
    table_inter(:,2) = (1:4^Ns)';
    table_inter(:,1) = M1*II_up';
    table_inter(:,3) = M1*II_down';
     
    [dummy,ind_s] = sort(table_inter(:,1));
     
    table = table_inter(ind_s,:);
     
    clear table_inter;
     
    table(:,4) = x11(ind_s11);
    table(:,5) = M_II_up(ind_s1);
     
%toc;
end
