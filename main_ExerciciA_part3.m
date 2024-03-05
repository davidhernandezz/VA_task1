clear; clc
%% Inputs
load('fe_model.mat')
dir_grvty = 1;          % 1, 2, 3 if the direction of gravity is X, Y or Z
write_h5file = false;

fixnodes = [10735 1 0;
            10735 2 0;
            10735 3 0;
            10735 4 0;
            10735 5 0;
            10735 6 0;
            
            13699 1 0;
            13699 2 0;
            13699 3 0;
            13699 4 0;
            13699 5 0;
            13699 6 0;
            
            16620 1 0;
            16620 2 0;
            16620 3 0;
            16620 4 0;
            16620 5 0;
            16620 6 0;
            
            19625 1 0;
            19625 2 0;
            19625 3 0;
            19625 4 0;
            19625 5 0;
            19625 6 0;
            
            22511 1 0;
            22511 2 0;
            22511 3 0;
            22511 4 0;
            22511 5 0;
            22511 6 0;
            
            4747  1 0;
            4747  2 0;
            4747  3 0;
            4747  4 0;
            4747  5 0;
            4747  6 0];

% Compute total number of nodes and DOFs
n_dof_total = size(K,1);
n_nodes = n_dof_total/6;

% Compute mass of each node
m = zeros(n_nodes,1);       % m: array with as many entries as nodes.
                                %    each entry corresponds to the mass of
                                %    that node.
for i = 1:n_nodes
    m(i,1) = M(6*(i-1)+1, 6*(i-1)+1);
end
m = m.*1000;        % mass is in [ton] in the original file -> [kg]

%% Boundary conditions
% Dirichlet
inD = zeros(size(fixnodes,1),1);    % array containing the DOFs in wich Diriclet conditions are applied.
uD = fixnodes(:,3);                 % the third column of fixnodes is uD

for i = 1:size(fixnodes,1)
    N = fixnodes(i,1);    
    inD(i) = 6*(N-1) + fixnodes(i,2);
end

a = 0;
inN = zeros(n_dof_total - size(fixnodes,1),1);
for i = 1:n_dof_total
    if ~ismember(i,inD)
        a = a + 1;
        inN(a) = i;
    end
end

% Neumann
F = zeros(n_dof_total,1);       % Array with as many entries as total num. of DOFs

for i = 1:size(inN,1)                % loop thtough inN (DOFs in wich Neumann conditions are applied)
    entry = inN(i);
    apply_force = false;        % a priori we assume there is no force in that DOF
    
    if mod(entry,6)-1 == 0            % if we are in a DOF multiple of 6, then:
        node = floor(entry/6)+1;        % the number that multiplied by 6 gives the DOF corresponds to the node of the DOF
        apply_force = true;
    end

    if apply_force                  % if we are in a DOF in wich no force should be applied, go to next
        if dir_grvty == 1
            dof = 6*(node-1) + 1;   % if gravity is in X-direction, DOFs where there is force: 1,7,13,etc
        elseif dir_grvty == 2
            dof = 6*(node-1) + 2;   % gravity in Y: DOFs where it actuates are 2,8,14, etc
        elseif dir_grvty == 3
            dof = 6*(node-1) + 3;   % grvty in Z, DOFs: 3,9,15, etc
        end
        node_mass = m(node);
        F(dof,1) = node_mass*9.81;  % [N] = [kg]*[m/sec^2]
    end
end

%% Partition
KDD = K(inD, inD);  % KRR
KDN = K(inD, inN);  % 
KNN = K(inN, inN);  % KLL
KND = K(inN, inD);
MNN = M(inN, inN);

FextN = F(inN,1);           % forces at Neumann nodes
FextD = F(inD,1);           % forces at Dirichlet nodes

selec=20;
tic
[eigvectors_restricted,eigenmodes_restricted] = eigs(KNN, MNN,selec,'smallestabs');
[eigvectors_free,eigenmodes_free] = eigs(K,M,selec,'smallestabs');
toc 

resonance_restricted_1=(eigenmodes_restricted).^(1/2);
resonance_restricted=(resonance_restricted_1)./(2*pi);

resonance_free_1=(eigenmodes_free).^(1/2);
resonance_free=(resonance_free_1)./(2*pi);
% Resultados: en efecto , los libres son 10^-7 veces mas pequeños que los
% restricted; no obstante los elementos restricted en la solución free no coinciden en la
% solución con los de restricted.


%% Ejercicio 4: Displacement frequency response function of the reference node; DUE TO UNIT ACCELERATION (1G) IN X DIRECTION FROM 0 TO 2000Hz (we will use 4000Hz as a computational limit) with 2% damping ratio
% We start by setting a frequency, going from 0Hz to 4000Hz in steps of 100
% -> We then solve the displacement using direct freq response with a new
% function only depending on the frequency

%We will also set the damping to 2%

phi=eigvectors_free;
phi_T = transpose(eigvectors_free);

mass_modal_B = transpose(phi)*M*phi;
stif_modal = transpose(phi)*K*phi; %Still to be multplied by omega
F_modal = transpose(phi)*F;
b = 0.02; %Critical Damping Factor

% mass_modal_A = phi_T*MNN;
% mass_modal_B = mass_modal_A*phi; %Still to be multplied by omega
% stif_modal_A = phi_T*KNN;
% stif_modal = mass_modal_A * phi; %Still to be multplied by omega
% F_modal = phi_T*FextN;
% b = 0.2; %Critical Damping Factor

counter = 1;
size_f_1=size(F);
size=size_f_1(1,1);
modal_analysis_matrix = zeros(size,selec);
omega_matrix = zeros(selec,1);

for omega=0:2:10
    omega_rad=omega*2*pi;
    mass_modal = mass_modal_B.*((-1)*omega_rad^2);
    Q_matrix = mass_modal + stif_modal;
    crit_damp = 1i * b * omega_rad;
    
    % Now length will be evaluated to extract those elements of the
    % diagonal to extract the value of eta_i and add the damping effect

    ref_A = length(Q_matrix);
    ref_B = ref_A(1,1);

    eta_vector = zeros(selec,1);

    for j=1:ref_B
        eta_vector(j,1) = F_modal(j,1)/(Q_matrix(j,j) + crit_damp);
    end
    
    % return of the eta_vector to cartesian coordiantes
    X = phi*eta_vector;
    modal_analysis_matrix(:,counter)=X;
    omega_matrix(counter,1) = omega;
    counter=counter+1;
end

%% Postprocessing of the displacements (storaged in modal_analysis_matrix)

% The objective here will be to extract the amplitude of the deformation of 
% the central node; for this the position in the modal_analysis_matrix of
% this node will be selected to then storage it into another position.

%position_searched = ((1304)*6)+1; 
position_searched = ((1305)*6); 
def = 0; %We will then compute the movement in 3D 
% of the central node with the components of X, Y and Z displacement.

for i=1:a
    if inN(i)==position_searched
        def = i;
        break
    else
    end
end

o = semilogx(omega_matrix(:,1), (modal_analysis_matrix(position_searched,:)), 'r');
grid on
o.LineWidth = 1.5;
title('x-Displacement Frequency Response Function (x-DFRF) of the Reference Node (Non-Restricted)')
xlabel('Frequency (Hz)')
xlim([0 2000]);
ylabel('Displacements (mm)')

