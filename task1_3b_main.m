clear
close all
clc


%% input data definition
% load K & M matrix
load('fe_model.mat');

% gravity acceleration vector
g = [0; 0; 0; 0; 0; 0];

% define n_dof 
dofs = 6;

% define the nodes where supports are located
n_supports = [];


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

d = eigs(K_nn,M(in_n, in_n),10,'smallestabs');

d = (sqrt(d)/(2*pi));