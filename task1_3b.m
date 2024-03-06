clear
close all
clc

%% DATA INPUT
% Load K & M matrix
load('fe_model.mat');

% Gravity acceleration vector
g = [0; 0; 0; 0; 0; 0];

% Define number of degrees of freedom at each node 
dofs = 6;

% Define the nodes where supports are located
% No supports are present, therefore the whole K and M matrix are Neumann 
n_supports = [];


%% COMPUTATIONS

% fix_nodes matrix

fix_nod = zeros(dofs * length(n_supports),3);

for i = 1:length(n_supports)
    fix_nod((dofs*i-(dofs-1):dofs*i),1) = n_supports(i);
    fix_nod((dofs*i-(dofs-1):dofs*i),2) = 1:dofs;
    fix_nod((dofs*i-(dofs-1):dofs*i),3) = 0;
end

% Dirichelt index vector
in_d = (fix_nod(:, 1) - 1) * dofs + fix_nod(:, 2);

% Dirichelt displacements vetor
u_d = fix_nod(:, 3);

% Neumann index matrix
in_n = setdiff(transpose(1:length(K)), in_d);

computedEigenmodes = 5; %5

% Eigenvectors and eigenvalues calculation
[V, D] = eigs(K,M,computedEigenmodes,'smallestabs');

% Frequency calculation
lambda = diag(D);
w = sqrt(lambda);
freq = (w/(2*pi));

eigModes(in_n,:) = V;

u = zeros(length(K)/dofs,dofs,computedEigenmodes);

% Displacements reshape for better visualitation
for i = 1:computedEigenmodes
    u(:,:,i) = transpose(reshape(eigModes(:,i), [dofs, length(K)/dofs]));
end

fillhdf('template.h5','1eigUncons.h5',u(:,:,1))
fillhdf('template.h5','5eigUncons.h5',u(:,:,5))