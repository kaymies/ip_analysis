function [M, C, J_c] = sub_Get_Com(m_i, c_i, j_i)
    % get COM 
    % in: m_i = [m_1, m_2, m_3, ...]
    %     c_i = [cx_1, cx_2, cx_3,...;
    %            cz_1, cz_2, cz_3, ...];
    if size(c_i,1) ~= 2
        if size(c_i,2) ==2
            c_i = c_i';
        else 
            error('c_i invalid for planar')
        end
    end

    nm = length(m_i);
    nc = size(c_i,2);

    if nm~=nc
        error('m_i, c_i not matching')
    end


    M = sum(m_i);
    if M == 0
        C = [sym(0);sym(0)];
    else
        C(1) = sum(m_i.*c_i(1,:))/M;
        C(2) = sum(m_i.*c_i(2,:))/M;
        C = C';
    end
    
    
    if nargin >2 && nargout >2
        J_c = 0;
        for i = 1:length(m_i)
            J_c = J_c + j_i(i)+m_i(i)*norm(C-c_i(:,i))^2;
        end
    end
    
end