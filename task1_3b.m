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

% Eigenvectors and eigenvalues calculation
[v, D] = eigs(K,M,11,'smallestabs');

% Frequency calculation
lamda = diag(D);
w = sqrt(lamda);
freq = (w/(2*pi));