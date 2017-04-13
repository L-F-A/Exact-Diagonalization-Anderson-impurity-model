function hij = Mat_element_AIM(bra,ket,Ns,ed,U,ee,V)
    if bra(length(bra)) ~= ket(length(ket))
        hij = 0;
    else
        
        Dket = ket(1)*ket(2);        
        nd = ket(1)+ket(2);
        ndup = ket(1);
        nddown = ket(2);
         ncup = zeros(1,(length(ket)-3)/2);
%         ncdown = zeros(1,(length(ket)-3)/2);
        nc = zeros(1,(length(ket)-3)/2);
        for r =1:length(nc)
            nc(r) = ket(2*r+1) + ket(2*r+2);
             ncup(r) = ket(2*r+1);
%             ncdown(r) = ket(2*r+2);
        end
        
%         ncup_tot = sum(ncup);
%         ncdown_tot = sum(ncdown);
        
        ele_test = sum((bra(1:length(bra)-1)-ket(1:length(ket)-1)).^2);
        ele = (ele_test == 0);
        
        for r=2:Ns
            
            if r == 2
                 pnci = 0;
            else
                pnci = sum(nc(1:r-2));
            end
           
           if  ket(1)*mod(ket(2*r-1)+1,2) == 0
                eleup_cd(r-1) = 0;
           else
                ket_up_cd = ket;
                %p = 2;
                p = nddown + pnci;

%               if r == 2
%                 p = 0;
%               elseif r == 3
%                 p = ncup(1);
%               else
%                 p = sum(ncup(1:r-2));
%               end
            
                ket_up_cd(1) = ket_up_cd(1)-1;
                ket_up_cd(2*r-1) = ket_up_cd(2*r-1)+1;
                eleup_cd(r-1) = ((-1)^p)*ket(1)*mod(ket(2*r-1)+1,2)*(sum((bra(1:length(bra)-1)-ket_up_cd(1:length(ket)-1)).^2) == 0);
           end
           
           if ket(2)*mod(ket(2*r)+1,2) == 0
                eledown_cd(r-1) = 0;
           else
                ket_down_cd = ket;
                %p = 2;
                p = 2*ndup + pnci + ncup(r-1);
            
%               if r == 2
%                 p = 2*sum(ncup);
%               elseif r == 3
%                 p = 2*sum(ncup) + ncdown(1);
%               else
%                 p = 2*sum(ncup) + sum(ncdown(1:r-2));
%               end


                ket_down_cd(2) = ket_down_cd(2)-1;
                ket_down_cd(2*r) = ket_down_cd(2*r)+1;
                eledown_cd(r-1) = ((-1)^p)*ket(2)*mod(ket(2*r)+1,2)*(sum( ( bra(1:length(bra)-1)-ket_down_cd(1:length(ket)-1)  ).^2 ) == 0);
           end
            
           if ket(2*r-1)*mod(ket(1)+1,2) == 0
                eleup_dc(r-1) = 0;
           else
                ket_up_dc = ket;
                %p = 2;
                p = nd + pnci;

%             if r == 2
%                 p = ket(1);
%             elseif r == 3
%                 p = ket(1)+ncup(1);
%             else
%                 p = ket(1) + sum(ncup(1:r-2));
%             end

                ket_up_dc(2*r-1) = ket_up_dc(2*r-1)-1;
                ket_up_dc(1) = ket_up_dc(1)+1;
                eleup_dc(r-1) = ((-1)^p)*ket(2*r-1)*mod(ket(1)+1,2)*(sum( ( bra(1:length(bra)-1)-ket_up_dc(1:length(ket)-1)  ).^2 ) == 0);
           end
           
           if ket(2*r)*mod(ket(2)+1,2) == 0
                eledown_dc(r-1) = 0;
           else
                ket_down_dc = ket;
                %p = 2;
                
                p = nd + pnci + ncup(r-1)+ndup;

%             if r == 2
%                 p = 2*sum(ncup)+ket(2);
%             elseif r == 3
%                 p = 2*sum(ncup)+ket(2)+ncdown(1);
%             else
%                 p = 2*sum(ncup)+ket(2) + sum(ncdown(1:r-2));
%             end

                ket_down_dc(2*r) = ket_down_dc(2*r)-1;
                ket_down_dc(2) = ket_down_dc(2)+1;
                eledown_dc(r-1) = ((-1)^p)*ket(2*r)*mod(ket(2)+1,2)*(sum( ( bra(1:length(bra)-1)-ket_down_dc(1:length(ket)-1) ).^2 ) == 0);
           end
            
        end
            hij = ( ed*nd + U*Dket + sum(ee.*nc) )*ele + sum(V.*( eleup_cd + eledown_cd + eleup_dc + eledown_dc ));
    end
    