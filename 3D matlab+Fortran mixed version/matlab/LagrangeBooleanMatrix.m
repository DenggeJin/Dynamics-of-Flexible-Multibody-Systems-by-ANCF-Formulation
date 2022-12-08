function Beq = LagrangeBooleanMatrix(dof_sys,elmnum_1,elmnum_2)


Beq = zeros(24,dof_sys,elmnum_1+elmnum_2);
for i = 1:elmnum_1
    for j = 1:24
        Beq(j,12*(i-1)+j,i) = 1;
    end
end
for i = elmnum_1+1:elmnum_1+elmnum_2
    for j = 1:24
        Beq(j,12*i+j,i) = 1;
    end
end

end