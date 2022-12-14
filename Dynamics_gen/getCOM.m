function [com_x, com_z, m_tot] = getCOM(list_mass, list_com_x, list_com_z)
    if((length(list_mass) ~= length(list_com_x))||(length(list_mass)~=length(list_com_z)))
        error('number of mass and number of length information do not match');
    end
    if(size(list_mass,1)~=1)
        error('put masses and coms in column by column');
    end
    
    n = length(list_mass);
    m_tot = sum(list_mass);
%     list_com_x = list_com(:,1);
%     list_com_z = list_com(:,2);

    com_x = list_mass'*list_com_x(:,1)/m_tot;
    com_z = list_mass'*list_com_x(:,2)/m_tot;
    
%     com = [com_x com_z]
    
end