% Compute edcosts with Amatrix
AF                                      = C_t{n};
AF(1+MMu:MMe,1+MMu:MMe)                 = sparse(I,I);
AF(MMM+1+MMu:MMM+MMe,MMM+1+MMu:MMM+MMe) = sparse(I,I);
AF(1+MMu:MMe,1+MMu:MMe)                 = spdiags((param.lx3_vec(1)+parent_ab_stay(1,1))*laction{n}(1:I,1),0,I,I)*Omega_t{n}; 
AF(MMM+1+MMu:MMM+MMe,MMM+1+MMu:MMM+MMe) = spdiags((param.lx3_vec(2)+parent_ab_stay(1,2))*laction{n}(1:I,2),0,I,I)*Omega_t{n};
AFT                  = AF';
AFT(1:MMu,:)         = 0;
AFT(MMe+1:MMM+MMu,:) = 0;
AFT(MMM+MMe+1:end,:) = 0;
if n<2
    edinflow_t(n,1)      = sum(AFT*(gg0.*bayyydelta(:)));
else
    edinflow_t(n,1)      = sum(AFT*(gold.*bayyydelta(:)));
end