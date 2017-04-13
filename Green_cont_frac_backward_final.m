function G = Green_cont_frac_backward_final(an,bn2,dd,EGS,direction,z)

    z = z+direction*EGS;

    P = 0;
    
    for r = 0:length(an)-1
        
        if r == length(an)-1
            P = ( z - direction*an(length(an)-r) ) - P;
            G = dd./P;
        else
            P = ( z - direction*an(length(an)-r) ) - P;
            P = bn2(length(bn2)-r)./P;
        end
    end

