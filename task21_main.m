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

y_displ_support = 1;
displ_supp = 2;

%% CALCULATIONS

fix_nod = fixnodes_p2_1(n_supports, dofs, displ_supp, y_displ_support);

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * 6 + fix_nod(:, 2);
% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

A = transpose(1:size(K, 1));
in_n = setdiff(A, in_d);

% gravity acceleration vector
g = [0; 9.81*10^3; 0; 0; 0; 0];

g_vect = repmat(g, 152340/6, 1);

F = M * g_vect;

F_n = F(in_n);

K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% displacements and forces vectors
u_n = K_nn\(F_n - K_nd * u_d);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/6]));
% displacements of the reference node
result = u(1305, :);
