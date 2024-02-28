function fix_nod = fixnodes_2a(n_supports, dofs, displ_supp, y_displ_support)
% fa la matriu de nodes, dof i desplaçament, pero imposa desplaçament en y
% per a un suport

fix_nod = zeros(dofs * length(n_supports),3);

for i = 1:length(n_supports)

    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;

end

fix_nod(dofs * displ_supp - 4, 3) = y_displ_support;