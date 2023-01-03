function u = structure(E, b, c, t, F_beam)

%-----------------------------------------------------------------------------
% Purpose:
%     Compute the structual analysis of a wing modeled as a rectangular beam
% 
% Synopsis:
%     u = structure(beam.E, beam.L, beam.b, beam.t, p)
% 
% Variable description (all in IU):
%   u: displacement vector for all the nodes
%   E: Elastic modulus
%   b: Beam's double length
%   c: Beam's width
%   t: Beam's thickness
%   p = vector of forces applied on the wing (deflection and twist)
%-----------------------------------------------------------------------------

E = E*1E9;
Ixx = c*t^3/12;                 % Beam's area moment of inertia
L = b/2;                        % Beams' length
ne = 100;                       % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof;                  
fdof = (3:DOF);                 % For a cantilever beam

% STIFFNESS MATRIX

coord_n = zeros(nn,2);                  % Nodal coordinates matrix
coord_n(:,1) = 0:L/ne:L;      % Only X coordinate is different than 0

connect_e = zeros(ne,2);                % Connectivity matrix of elements through nodes
connect_e(:,1) = 1:1:ne;
connect_e(:,2) = 2:1:nn;

K = zeros(DOF);                         % Initilization of the stiffness matrix

for e = 1:ne
    index = connect_e(e,:);                                    
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2 - x1;                                               % Length of the element

    dofe = [index(1)*dofn-1 index(1)*dofn...                    % DOF of each element
            index(2)*dofn-1 index(2)*dofn];

% STIFFNESS
     
    k = E*Ixx/Le^3;
    Kef = k*[12 6*Le -12 6*Le;...
             6*Le 4*Le^2 -6*Le 2*Le^2;...
             -12 -6*Le 12 -6*Le;...
             6*Le 2*Le^2 -6*Le 4*Le^2];

    K(dofe, dofe) = K(dofe, dofe) + Kef;

end

% Solve the static problem: [K]*{u} = {p}

K = K(fdof,fdof);
F_beampc = F_beam(fdof);

u = K\F_beampc;

end
