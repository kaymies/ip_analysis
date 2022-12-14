function inv_T = inv_Trans_mat(T)
    R = T(1:2,1:2);
    t = T(1:2, 3);
    
    inv_T(1:2, 1:2) = R';
    inv_T(1:2,3) = -R'*t;
    inv_T(3,1:2) = zeros(1,2);
    inv_T(3,3) = 1;
    simplify(inv_T);
end