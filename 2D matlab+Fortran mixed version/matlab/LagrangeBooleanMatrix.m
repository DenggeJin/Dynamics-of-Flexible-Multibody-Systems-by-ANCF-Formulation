function Beq = LagrangeBooleanMatrix(dof_sys,elmnum_1,elmnum_2)


Beq = zeros(12,dof_sys,elmnum_1+elmnum_2);
for i = 1:elmnum_1
    for j = 1:12
        Beq(j,6*(i-1)+j,i) = 1;
    end
end
for i = elmnum_1+1:elmnum_1+elmnum_2
    for j = 1:12
        Beq(j,6*i+j,i) = 1;
    end
end

end