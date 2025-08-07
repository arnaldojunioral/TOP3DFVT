% TOPOLOGY OPTIMIZATION OF 3D ELASTIC STRUCTURES BY ZEROTH ORDER
% FINITE-VOLUME THEORY

function TOP3DFVT(n1, n2, n3, volfrac)
%% ____________________________________________________________USER-DEFINED
[L, H, B] = deal(n1, n2, n3);                                             % Beam dimensions
[E0, Emin, nu]  = deal(1, 1e-9, 0.3);                                      % Material properties
P  = -1;                                                                   % Applied load
model = 'GSS';                                                             % Set to 'GSS' for grayscale suppression (linear interpolation)
eta = 1/3;                                                                 % Filtering damping factor
tol = 1e-3;                                                                % Convergence tolerance
maxit = 1000;                                                              % Maximum number of iterations
threshold = 0.01;                                                          % Density threshold for display
displayflag = true;                                                        % Display during optimization (true = live updates)
pb = 'Cantilever';                                                         % Problem type

%% ________________________________________________________________________

% Interpret optional inputs depending on model
if any(strcmpi(model, {'SIMP', 'RAMP'}))
    grayscale = false;

elseif strcmpi(model, 'GSS')
    grayscale = true;

else
    error('Unknown material model: %s', model);
end

%% _____________________________________ INTERPOLATION / PENALIZATION MODEL
switch upper(model)
    
    case 'SIMP'  % Solid Isotropic Material with Penalization
        MatInt = @(p, x) {Emin + (E0-Emin) * x.^p, ...
            p * (E0-Emin) * x.^(p-1)};

    case 'RAMP'  % Rational Approximation of Material Properties
        MatInt = @(p, x) {Emin + (E0-Emin)*x./(1 + p*(1 - x)), ...
            (1 + p)*(E0 - Emin)./(1 + p*(1 - x)).^2};

    case 'GSS'   % Grayscale suppression – linear interpolation (Voigt model)
        MatInt = @(~, x) {Emin + (E0 - Emin)*x, (E0 - Emin) * ones(size(x))};
        penal = 1; % dummy value just to enter the loop once

    otherwise
        error('Unknown material model: %s. Choose SIMP, RAMP, or GSS.', model);
end

%% ________________________________________________________DESIGN VARIABLES

% subvolume dimensions
[l, h, b] = deal(L/n1, H/n2, B/n3);

% Initialize uniform material density
x = volfrac * ones(n1, n2, n3);

% Define active design region depending on the problem type
switch pb
    
    case 'Bridge'
        % Fix solid subvolumes at the far end of y-direction
        x(:, end-4:end, :) = 1;

        % Create logical mask for active design variables
        mask = true(n1, n2, n3);
        mask(:, end-4:end, :) = false;

        % Extract active indices for optimization
        active = find(mask);

        % Regularize density to match volume constraint
        x = regularizeDensity(x, active, volfrac);

    otherwise
        % All subvolumes are active (default case)
        active = 1:numel(x);
end

% Volume constraint gradient
dvdx = repmat(l * h * b, [n1, n2, n3]);

%% ___________________________________PREPARE FINITE-VOLUME THEORY ANALYSIS

[dof, ndof, iK, jK] = DOFassembly(n1, n2, n3);                             % degrees of freedom
K0 = LocalStiffMatrix(E0, nu, l, h, b);                                    % local stiffness matrix

% BOUNDARY CONDITIONS: global force and fixed dofs (support)
switch pb

    case 'Cantilever'
        F = sparse(dof(n1:n1 * n2:end, 5)', 1,  P, ndof, 1);
        supp = unique(dof(1:n1:end, 1:3));

    case 'MBB'
        F    = sparse(dof(n1 * (n2 - 1) + 1:n1 * n2:end, 11)', 1,  P, ndof, 1);
        supp = unique([unique(dof(1:n1:end, [1,3])); dof(n1:n1 * n2:end, 8)]);
        % supp = unique([unique(dof(1:n1:end, [1,3])); dof(n1:n1 * n2 * (n3-1):end, 8)]);

    case 'Bridge'
        % Map subvolumes on the top layer (j = n2)
        [i, k] = ndgrid(1:n1, 1:n3);
        subvol = i + (n2 - 1) * n1 + (k - 1) * n1 * n2;

        % Load vector (column 11 of the 'dof' matrix)
        F = sparse(dof(subvol(:), 11)', 1,  P, ndof, 1);

        % Select base element indices along x (i) for supports at j = 1
        id = ceil(n1 / 4);
        base_ids = [id, n1 - id + 1];

        % Propagate base elements along the z-direction (k)
        layer = (0:n3-1)' * n1 * n2;
        subvol = base_ids + layer;

        % Fixed DOFs (columns 7 to 9 of 'dof' matrix)
        supp = unique(dof(subvol(:), 7:9));
end

% Global stiffness matrix assemblage
fdof = setdiff(dof(:), supp(:));                                           % free degrees of freedom
StiffnessAssemblage = @(sK) sparse(iK(:), jK(:), sK(:), ndof, ndof);       % global stiffness matrix

%% ____________________________________________________OPTIMIZATION PROCESS
t_start = tic;
for i = 1:length(penal(:))

    [change, loop, qmax, q, xPhys] = deal(1.0, 0, 1/eta, 1.0, x);
    fprintf('\nPenalty factor: %1.2f\n', penal(i));                        % print current penalty factor

    while (change > tol && loop < maxit), loop = loop+1;                   % start optmization process

        % GRAY SCALE CONTROLLER
        q = 1 + grayscale * (loop >= 20) * (min(qmax, 1.02 * q) - 1);

        % FINITE-VOLUME THEORY ANALYSIS
        Mat = MatInt(penal(i), xPhys);                                     % material interpolation
        E = Mat{1}; dEdx = Mat{2};
        sK = K0(:) * E(:)';                                                % stiffness interpolation
        K = StiffnessAssemblage(sK); K = (K+K')/2;                         % assemblage of the stiffness matrix
        U = Tikhonov(fdof, K, F);                                          % compute global displacements employing Tikhonov strategy

        % COMPLIANCE AND SENSITIVITY
        fe = reshape(sum((U(dof) * K0) .* U(dof), 2),[n1, n2, n3]);        %
        f  = sum(sum(sum(E .* fe)));                                       % objective function: structural compliance
        dfdx = -dEdx .* fe;                                                % sensitivity

        % UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2 - l1)/(l1 + l2) > 1e-3
            lmid = 0.5*(l2 + l1);
            beta = (-dfdx ./ dvdx / lmid).^eta;
            xnew = x;
            xnew(active) = max(0, max(x(active)-move, min(1, min(x(active)+move, (x(active).*beta(active)).^q))));
            xPhys(active) = xnew(active);
            if mean(xPhys(:)) > volfrac
                l1 = lmid;
            else
                l2 = lmid;
            end
        end

        % Update design variable
        change = max(abs(xnew(:) - x(:)));
        x = xnew;

        % PRINT RESULTS
        fprintf('It: %i\tObjec.: %1.2f\tVol.: %1.2f\tGSS: %1.2f\tChange: %1.3f\n',...
            loop, f, mean(xPhys(:)), q, change);

        % PLOT DENSITIES
        if displayflag, clf; display3D(xPhys, l, h, b, threshold); end

        % VOLUME CONTROLLER
        if abs(mean(xPhys(:)) - volfrac) > 0.001
            warning('Volume constraint deviated: current = %.4f, target = %.4f', ...
                mean(xPhys(:)), volfrac);
            return;
        end
    end
end

% BLACK-AND-WHITE FRACTION
bw = nnz(xPhys < 0.001+1e-6 | xPhys > 1.0-1e-6) / numel(xPhys);
fprintf('B&W: %1.3f\n', bw);

% CHECK SYMMETRY ACROSS XY-PLANE
symmetry_error = max(abs(xPhys - flip(xPhys, 3)), [], 'all');
is_symmetric_xy = symmetry_error < 1e-3;
fprintf('Symmetry across xy-plane: %s\n', string(is_symmetric_xy));

% PLOT DESIGN
if ~displayflag, clf; display3D(xPhys, l, h, b, threshold); end

% PROCESSING TIME
t_total = toc(t_start);
fprintf('Total processing time: %0.2f seconds \n', t_total);

%% _____________________________________________ORDENING DEGREES OF FREEDOM
function [dof, ndof, iK, jK] = DOFassembly(n1, n2, n3)

% Indices for subvolumes:
[i, j, k] = ndgrid(1:n1, 1:n2, 1:n3);

% Number of horizontal local_faces:
Nhf = n1 * n3 * (n2 + 1);

% Number of vertical local_faces in the x-direction:
Nvfx = (n1 + 1) * n3 * n2;

% Subvolume local_faces:
fc1 = Nhf + i(:) + (k(:) - 1) * (n1 + 1) + (j(:) - 1) * n3 * (n1 + 1);     % Left lateral local_faces
fc2 = Nhf + i(:) + 1 + (k(:) - 1) * (n1 + 1) + (j(:) - 1) * n3 * (n1 + 1); % Right lateral local_faces
fc3 = i(:) + (k(:) - 1) * n1 + (j(:) - 1) * n1 * n3;                       % Bottom local_faces
fc4 = i(:) + (k(:) - 1) * n1 + j(:) * n1 * n3;                             % Top local_faces
fc5 = Nhf + Nvfx + i(:) + (k(:) - 1) * n1 + (j(:) - 1) * n1 * (n3 + 1);    % Back local_faces
fc6 = Nhf + Nvfx + i(:) + k(:) * n1 + (j(:) - 1) * n1 * (n3 + 1);          % Front local_faces

% Final combination of local_faces:
local_faces = [fc1, fc2, fc3, fc4, fc5, fc6];

% Degrees of fereedom
dof = zeros(n1 * n2 * n3, 18);
dof(:, 1:3:16) = 3 * local_faces - 2;
dof(:, 2:3:17) = 3 * local_faces - 1;
dof(:, 3:3:18) = 3 * local_faces;

ndof = max(dof(:));                                                        % total number of degrees of freedom
iK = kron(dof, ones(18, 1))';                                              % iK - line mapping indice to the global stiffness location
jK = kron(dof, ones(1, 18))';                                              % jK - column mapping indice to the global stiffness location

%% ____________LOCAL STIFFNESS MATRIX FOR ZEROTH ORDER FINITE VOLUME THEORY
function K0 = LocalStiffMatrix(E0, nu, l, h, b)

% Constitutive matrix
factor = E0 / ((1 + nu) * (1 - 2 * nu));
C0 = factor * [
    1 - nu, nu,     nu,     0,         0,         0;
    nu,     1 - nu, nu,     0,         0,         0;
    nu,     nu,     1 - nu, 0,         0,         0;
    0,      0,      0,      (1 - 2*nu)/2, 0,         0;
    0,      0,      0,      0,         (1 - 2*nu)/2, 0;
    0,      0,      0,      0,         0,         (1 - 2*nu)/2];

% Normal vectors
N = zeros(18,36);
N1 = [-1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -1; 0, 0, 0, 0, -1, 0];
N3 = [0, 0, 0, 0, 0, -1; 0, -1, 0, 0, 0, 0; 0, 0, 0, -1, 0, 0];
N5 = [0, 0, 0, 0, -1, 0; 0, 0, 0, -1, 0, 0; 0, 0, -1, 0, 0, 0];

N(1:3, 1:6)     = N1;
N(4:6, 7:12)    = -N1;
N(7:9, 13:18)   = N3;
N(10:12, 19:24) = -N3;
N(13:15, 25:30) = N5;
N(16:18, 31:36) = -N5;

% Kinematic matrices:
a = repmat(eye(3), 6, 1);
A = sparse(18, 18);

A(1,1)  = -l/2;    A(1,4)  = l^2/4;
A(2,7)  = A(1,1);  A(2,10) = A(1,4);
A(3,13) = A(1,1);  A(3,16) = A(1,4);
A(4,1)  = -A(1,1); A(4,4)  = A(1,4);
A(5,7)  = -A(1,1); A(5,10) = A(1,4);
A(6,13) = -A(1,1); A(6,16) = A(1,4);

A(7,2)   = -h/2;    A(7,5)   = h^2/4;
A(8,8)   = A(7,2);  A(8,11)  = A(7,5);
A(9,14)  = A(7,2);  A(9,17)  = A(7,5);
A(10,2)  = -A(7,2); A(10,5)  = A(7,5);
A(11,8)  = -A(7,2); A(11,11) = A(7,5);
A(12,14) = -A(7,2); A(12,17) = A(7,5);

A(13,3)  = -b/2;     A(13,6)  = b^2/4;
A(14,9)  = A(13,3);  A(14,12) = A(13,6);
A(15,15) = A(13,3);  A(15,18) = A(13,6);
A(16,3)  = -A(13,3); A(16,6)  = A(13,6);
A(17,9)  = -A(13,3); A(17,12) = A(13,6);
A(18,15) = -A(13,3); A(18,18) = A(13,6);

% Strain/Displacement operator per local_faces: strain = E * Wi
KinematicRelation = @(x, y, z) [1, 0, 0, 3*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*y, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*z;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*z, 0, 1, 0, 0, 3*y, 0;
    0, 0, 1, 0, 0, 3*z, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*x, 0, 0;
    0, 1, 0, 0, 3*y, 0, 1, 0, 0, 3*x, 0, 0, 0, 0, 0, 0, 0, 0];

E1 = KinematicRelation(-l/2, 0, 0);
E2 = KinematicRelation(l/2, 0, 0);
E3 = KinematicRelation(0, -h/2, 0);
E4 = KinematicRelation(0, h/2, 0);
E5 = KinematicRelation(0, 0, -b/2);
E6 = KinematicRelation(0, 0, b/2);
E = [E1; E2; E3; E4; E5; E6];

% Constitutive matrix: 36 by 36 elements
C = zeros(36, 36);

% Use a loop to place C0 in the diagonal blocks
idx = 1:6:36;  % Indices for each block
for i = idx
    C(i:i+5, i:i+5) = C0;  % Assign C0 to each block
end

% Operator that relates the local surlocal_faces-averaged tractions x Wi
B = N * C * E;

% Eficcienty matrix inverse
identI = eye(18, 18);
invA = A \ identI;

% Summation components of matrix (B * S) per local_faces
sumB = (B(1:3, :) + B(4:6, :)) * b * h + ...
    (B(7:9, :) + B(10:12, :)) * l * b + ...
    (B(13:15, :) + B(16:18, :)) * l * h;

% Matrices a_barra and A_barra
ab = (sumB * invA * a) \ (sumB * invA);
Ab = invA * (identI - a * ab);

% Pseudo local stiffness matrix
K0 = B * Ab;

% Create an 18x18 matrix of zeros for the subvolume local_facess area
S = zeros(18, 18);

% Fill the first 6x6 block with b * h on the diagonal
S(1:6, 1:6) = b * h * eye(6);

% Fill the next 6x6 block with l * b on the diagonal
S(7:12, 7:12) = l * b * eye(6);

% Fill the last 6x6 block with l * h on the diagonal
S(13:18, 13:18) = l * h * eye(6);

% Modified local stiffness matrix
K0 = S * K0;

%% __________________________________________TIKHONOV REGULARIZATION SCHEME
function Unew = Tikhonov(fdof, K, F)

% Define tolerance
tol = 1e-6;

% Number of degrees of freedom
ndof = length(K);

% Initial regularization parameter
lambda0 = 10.^(-12);
lambda = lambda0 .* trace(K) / ndof;

% Initialize variables
U = zeros(ndof, 1); Unew = U;
error = 1; iter = 1;

while error > tol

    % Tikhonov regularization
    Knew = K + lambda * speye(ndof);

    % Solve the system of equations
    Unew(fdof) = Knew(fdof, fdof) \ F(fdof);

    % Compute the relative error
    error = norm(K(fdof, fdof) * Unew(fdof) - F(fdof)) / norm(F(fdof));

    % Update the displacement vector
    U(fdof) = Unew(fdof);

    % iteration
    iter = iter + 1;
end

%%_______________________________________________________REGULARIZE DENSITY
function x = regularizeDensity(x, active, volfrac)

% Scale active elements to enforce mean(x(:)) = volfrac
total_elements = numel(x);
fixed_sum = sum(x(:)) - sum(x(active));  % = 0 if all elements are active

% Compute the target sum for the active region to enforce global mean = volfrac
target_total = volfrac * total_elements;
target_active_sum = target_total - fixed_sum;

% Apply scaling factor to active region
x(active) = x(active) * target_active_sum / sum(x(active));

%% _________________________________________PLOTTING THE OPTIMIZED TOPOLOGY
function display3D(rho, l, h, b, threshold)

% Mesh discretization
[n1, n2, n3] = size(rho); nsv = n1 * n2 * n3;

% Predefined local faces
local_faces = [1 4 3 2; 1 2 6 5; 2 3 7 6; 3 4 8 7; 1 5 8 4; 5 6 7 8];

% Preallocate storage for faces, vertices and colors
Vertices = zeros(8 * nsv, 3);
Faces = zeros(6 * nsv, 4);
colors = zeros(6 * nsv, 3);

% Define base coordinates for the subvolume
vx = [0, 1, 1, 0, 0, 1, 1, 0];
vy = [0, 0, 1, 1, 0, 0, 1, 1];
vz = [0, 0, 0, 0, 1, 1, 1, 1];

% Loop over all voxels
idV = 0; idF = 0;
for k = 1:n3
    for j = 1:n2
        for i = 1:n1

            % Check if the density of the subvolume is greater than the threshold
            if (rho(i, j, k) > threshold)

                % Compute subvolume vertices
                x = (i - 1) * l + vx * l;
                y = (j - 1) * h + vy * h;
                z = (k - 1) * b + vz * b;

                coord = [x(:), y(:), z(:)];

                % Color based on density (lighter for lower density)
                val = 0.4 + 0.6 * (1 - rho(i, j, k));

                % Add this voxel's data
                Vertices(idV + (1:8), :) = coord;
                Faces(idF + (1:6), :) = local_faces + idV;
                colors(idF + (1:6), :) = repmat([val val val], 6, 1);

                % Update indices
                idV = idV + 8;
                idF = idF + 6;
            end
        end
    end
end

% Removing subvolumes with densities lower than the predefined threshold
Vertices = Vertices(1:idV, :);
Faces = Faces(1:idF, :);
colors = colors(1:idF, :);

% Swap columns [2 3]
Vertices(:, [3, 2]) = Vertices(:, [2, 3]);

% Plot all patches at once
patch('Faces', Faces, 'Vertices', Vertices, 'FaceVertexCData', colors, ...
    'FaceColor', 'flat');
axis equal; axis tight; axis off; box on;
view([30, 20]);
set(gcf, 'Color', 'w');

% Configure figure window
fig = figure(1);
set(fig, 'Name', 'Top3DFVT', 'NumberTitle', 'off');

%% _____________________________________________________________________END