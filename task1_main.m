clear
close all
clc


%% input data definition
% load K & M matrix
load('fe_model.mat');

% gravity acceleration vector
g = [0; 9.81*10^3; 0; 0; 0; 0];

% define n_dof 
dofs = 6;

% define the nodes where supports are located
n_supports = [10735; 13699; 16620; 19625; 22511; 4747];


%% Preliminary computations
% fix_nodes matrix
fix_nod = fixnodes(n_supports, dofs);

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);
% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

in_n = setdiff(transpose(1:length(K)), in_d);

g_vect = repmat(g, length(M)/dofs, 1);

Fext = M * g_vect;

F_n_ext = Fext(in_n);
F_d_ext = Fext(in_d);

K_nn = K(in_n, in_n);
K_dd = K(in_d, in_d);
K_nd = K(in_n, in_d);
K_dn = K(in_d, in_n);

% displacements and forces vectors
u_n = K_nn\(F_n_ext - K_nd * u_d);

u(in_n, 1) = u_n;
u(in_d, 1) = u_d;

u = transpose(reshape(u, [dofs, length(K)/dofs]));

F_d = K_dd*u_d + K_dn*u_n;

F(in_n) = F_n_ext;
F(in_d) = F_d + F_d_ext;
F = transpose(reshape(F, [dofs, length(K)/dofs]));

%fillhdf('template.h5','example.h5',u)

sumReactions(1) = sum(F(:,1));
sumReactions(2) = sum(F(:,2));
sumReactions(3) = sum(F(:,3));
sumReactions(4) = sum(F(:,4));
sumReactions(5) = sum(F(:,5));
sumReactions(6) = sum(F(:,6));

F_d = transpose(reshape(F_d, [dofs, length(F_d)/dofs]));
% F_n_ext = transpose(reshape(F_n_ext, [dofs, length(F_n_ext)/dofs]));
% F_d_ext = transpose(reshape(F_d_ext, [dofs, length(F_d_ext)/dofs]));
% 
% sumreactions = sum(F_d(:,2));
% sumreactionsM = sum(F_n_ext(:,2)) + sum(F_d_ext(:,2));
% sumreactionsM = sum(F_n_ext(:,2)) + sum(F_d_ext(:,2));
% 
% result = sumreactionsM + sumreactions

% postprocessar resultats
    % obrir meta
    % afegir arxiu h5 al meta
    % plotejar desplaçaments --> pestaña results, read 
    % fringe --> activar colors de tensions del model
    % armonic 



