function coe_el = mee2coe(mee_el)
    p = mee_el(1);
    f = mee_el(2);
    g = mee_el(3);
    h = mee_el(4);
    k = mee_el(5);
    L = mee_el(6);
        
    e = sqrt(f^2+g^2);
    a = p/(1-e^2);
    
    if f == 0 && g == 0 || e ==0
        w_ohm = 0;
    else
        w_ohm = atan2(g/e,f/e);
    end
        
    if h == 0 && k == 0
        i = 0;
        ohm = 0;
    else
        tani_2 = sqrt(h^2+k^2);
        i = atan(tani_2)*2;
        ohm = atan2(k/(tan(i/2)),h/(tan(i/2)));
    end

    

    w = w_ohm-ohm;
    nu = L - (w_ohm);
    
    coe_el = [a,e,i,w,ohm,nu];
end