clear
close all
clc

%% DATA INPUT
% load K & M matrix
load('fe_model.mat');

% define n_dof
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];
% define the referecne node
ref_node = 1305;

y_displ_support = 1; %unit displacement

result = zeros(length(n_supports),dofs);

for displ_supp = 1:length(n_supports) %computes the reference point displacement for a given displacement in the y axis of each of the supports

    %% CALCULATIONS

    fix_nod = fixnodes_2a(n_supports, dofs, displ_supp, y_displ_support);

    % Dirichelt index vector
    in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
    % Dirichelt displacements vetor
    u_d = fix_nod(:, 3);

    in_n = setdiff(transpose(1:length(K)), in_d);

    % gravity acceleration vector
    g = [0; 9.81*10^3; 0; 0; 0; 0];

    g_vect = repmat(g, length(K)/dofs, 1);

    F = M * g_vect;

    F_n = F(in_n);

    K_nn = K(in_n, in_n);
    K_dd = K(in_d, in_d);
    K_nd = K(in_n, in_d);
    K_dn = K(in_d, in_n);

    % displacements and forces vectors
    u_n = K_nn\(F_n - K_nd * u_d);

    u = zeros(length(K),1);

    u(in_n, 1) = u_n;
    u(in_d, 1) = u_d;

    u = transpose(reshape(u, [dofs, length(K)/dofs]));
    % displacements of the reference node
    result(displ_supp,:) = u(ref_node, :);

end
