function [m_tot, length_tot, com_tot, j_tot] = sub_paramconversion(list_mass, list_length, list_com, list_gyration, ind_flip)
    
    n = length(list_mass);
    m_tot = sum(list_mass);
    
    
    
    
    
    tmp = 0;
    com_length = 0;
    list_com_cum = zeros(n,1);
    for i = 1:n        
        list_com_cum(i) = com_length + list_com(i);
        tmp = tmp + list_mass(i)*(list_com_cum(i));
        com_length = com_length + list_length(i);        
    end
    com_tot = tmp/m_tot;
    length_tot = com_length;
    
    
    
    
    tmp = 0;
    for i = 1:n
        tmp = tmp + list_mass(i)*( ...
            list_gyration(i)^2 ...
            +(com_tot-list_com_cum(i))^2 ...
        );              
    end
        
    j_tot = tmp; % about segment total COM
    
    if ind_flip
        com_tot = length_tot - com_tot; % to flip the original coordinate to my coordinate
    end
    
end


