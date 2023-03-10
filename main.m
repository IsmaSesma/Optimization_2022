%--------------------- MAIN - TOCDA-----------------------               
%                         1 - 2023  
%
% Celia Ayuso Ramos
% Jaime Gómez Benavente
% Ángel Ladrón Córdoba
% Hugo López Fernandez
% Ismael Rodríguez Sesma

clc; clear; close all
warning('off','all')

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

%% Initial values

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")
    
M = 2; % Mach de vuelo para crucero
M_t = 1.25; % Mach de vuelo para giro
AOA = 2; % Ángulo de ataque para crucero [deg]
AOA_t = AOA; % Ángulo de ataque para giro [deg]
Sweep = 50; % Flecha [deg]
Wingspan = 12; % Envergadura [m]
hft = 30000; % Altitud para crucero [ft]
hft_t = 3300; % Altitud para giro [ft]
Fm = 8000; % Masa combustible [kg]
E = 75E9;   % Modulo de Young

%% Restricciones

% Restricciones de desigualdad

A = [1 0 0 0 0 0 0 0 0 0; ...
     -1 0 0 0 0 0 0 0 0 0; ...
     0 1 0 0 0 0 0 0 0 0; ...
     0 -1 0 0 0 0 0 0 0 0; ...
     0 0 0 0 1 0 0 0 0 0; ...
     0 0 0 0 -1 0 0 0 0 0; ...
     0 0 0 0 0 1 0 0 0 0; ...
     0 0 0 0 0 -1 0 0 0 0; ...
     0 0 0 0 0 0 0 0 0 1; ...
     0 0 0 0 0 0 0 0 0 -1];

b = [2.5 -1.5 1.5 -1 70 -50 25 -10 135 -70];

% Restricciones de igualdad

Aeq = [0 0 1 0 0 0 0 0 0 0; ...
       0 0 0 1 0 0 0 0 0 0 ; ...
       0 0 0 0 0 0 1 0 0 0; ...
       0 0 0 0 0 0 0 1 0 0; ...
       0 0 0 0 0 0 0 0 1 0];

beq = [AOA AOA_t hft hft_t Fm/1E3];

% Valores iniciales externos

vect_init_mono = [M, M_t, AOA, AOA_t, Sweep, Wingspan, hft, hft_t, Fm/1E3, E/1E9];

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

%% Simulador test
 
% Aerodinámica
[CL, ~, CD, Lift, Drag, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
t = 0.07*c;     % Espesor del ala

[CL_t, ~, CD_t, Lift_t, Drag_t, F_beam_t, c_t, S_t, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
t_t = 0.07*c_t;     % Espesor del ala

% Estructura
u = structure(E, Wingspan, c, t, F_beam);
u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);

% Actuaciones y propulsión
[R, Aut, To, Rt, To_t, gforce, mu_deg, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm, v_t, Lift_t, Drag_t, Lift, valores_T,...
                                                                  valores_SFC, valores_T_af, valores_SFC_af);

% Costes
[coste_motor, precio_mat, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, To_min, Fm, valores_E,...
                                                                       precio_materiales, E/1E9, S, precio_combustible);


%% Optimización mono-objetivo

% Gradiente
options = optimoptions('fmincon','Display','iter','MaxIterations',6, 'Algorithm','sqp');
[X, F, exitflag, output, lambda, grad, hessian] = fmincon(@coste_monoobjetivo, vect_init_mono, A, b', Aeq, beq', [], [], [], options);

% Análisis de sensibilidad
changeable = [1,2,5,6,10];

for i=1:length(changeable)
    upperb = X(changeable(i))*1.05;
    lowerb = X(changeable(i))*0.95;
    
    sens_X_vec = linspace(lowerb, upperb, 50);
    values = zeros(50,1);
    for j =1:50
        input_vec_mod = X;
        input_vec_mod(changeable(i)) = sens_X_vec(j);
        
        values(j) = coste_monoobjetivo(input_vec_mod);
    end
    
    figure(i)
    hold on
    plot(sens_X_vec, values)
end

%% Optimización heurística

% Reescalado

rescale_init(:,1) = vect_init_mono(:,1)*1E3;
rescale_init(:,2) = vect_init_mono(:,2);
rescale_init(:,3) = vect_init_mono(:,3)*1E2;
rescale_init(:,4) = vect_init_mono(:,4);
rescale_init(:,5) = vect_init_mono(:,5)*1E2;
rescale_init(:,6) = vect_init_mono(:,6)*sqrt(1E5);
rescale_init(:,7) = vect_init_mono(:,7);
rescale_init(:,8) = vect_init_mono(:,8);
rescale_init(:,9) = vect_init_mono(:,9);
rescale_init(:,10) = vect_init_mono(:,10);

% Restricciones de desigualdad

b_heu = [2.5*1E3 -1.5*1E3 1.5 -1 70*1E2 -50*1E2 25*sqrt(1E5) -10*sqrt(1E5) 135 -70];

% Restricciones de igualdad

beq_heu = [AOA*1E2 AOA_t hft hft_t Fm/1E3];

% Optimiación con reescalado -- gradiente

[X_heu, F_heu, exitflag_heu, output_heu, lambda_heu, grad_heu, hessian_heu] = fmincon(@coste_heu, rescale_init, A, b_heu', Aeq, beq_heu', [], [], [], options);

% Optimización -- heurística

    % Genetic Algorithm
options_ga = optimoptions('ga','Display','iter');
[X_ga, F_ga, exitflag_ga, output_ga, population_ga, scores_ga] = ga(@coste_heu, 10, A, b_heu, Aeq, beq_heu, [], [], [], [], options_ga);

 

    % Genetic Algorithm Modificado
options_gam = optimoptions('ga','Display','iter','MaxGenerations',5,'HybridFcn','fmincon','PopulationSize',20);
[X_gam, F_gam, exitflag_gam, output_gam, population_gam, scores_gam] = ga(@coste_heu, 10, A, b_heu, Aeq, beq_heu, [], [], [], [], options_ga);

 

    % Pattern Search
options_PS = optimoptions("patternsearch", 'Display','iter');
[X_PS, F_PS, exitflag_PS, output_PS] = patternsearch(@coste_heu, rescale_init, A, b_heu, Aeq, beq_heu, [], [], [], options_PS);

 

    % Pattern Search Modificado
options_PSm = optimoptions("patternsearch", 'Display','iter','MaxIterations',100,'MeshContractionFactor',0.2,'MeshExpansionFactor',2.5);
[X_PSm, F_PSm, exitflag_PSm, output_PSm] = patternsearch(@coste_heu, rescale_init, A, b_heu, Aeq, beq_heu, [], [], [], options_PS);

 

    % Particle Swarm
options_PSW = optimoptions("particleswarm", 'Display','iter');
[X_PSW, F_PSW, exitflag_PSW, output_PSW] = particleswarm(@coste_heu, 10, [1.5*1E3;1;3*2.52856e+04;2;50*1E2;10*sqrt(1E5);30000;3300;8;70],...
                                            [2.5*1E3;1.5;3*2.52856e+04;2;70*1E2;25*sqrt(1E5);30000;3300;8;135], options_PSW);

    % Particle Swarm Modificado
options_PSWm = optimoptions("particleswarm", 'Display','iter','MaxIterations',10,'HybridFcn','fmincon','SwarmSize',20);
[X_PSWm, F_PSWm, exitflag_PSWm, output_PSWm] = particleswarm(@coste_heu, 10, [1.5*1E3;1;3*2.52856e+04;2;50*1E2;10*sqrt(1E5);30000;3300;8;70],...
                                            [2.5*1E3;1.5;3*2.52856e+04;2;70*1E2;25*sqrt(1E5);30000;3300;8;135], options_PSW);

%% Optimización Multiobjetivo

    % Genetic Algorithm Mutliobj
options_mo = optimoptions('gamultiobj', 'Display', 'iter', 'PlotFcn',{'gaplotpareto'}, 'PopulationSize', 30, 'MaxGenerations', 8, 'ParetoFraction', 0.7);
[x_gamo,fval_gamo,exitflag_gamo,output_gamo,residuals_gamo] = gamultiobj(@coste_mo, 10, A, b_heu, Aeq, beq_heu, [], [], options_mo);


    % paretosearch
options_mo = optimoptions('paretosearch', 'Display', 'iter', 'PlotFcn',{'psplotparetof'}, 'ParetoSetSize', 20, 'MaxIterations', 10);
[x_p,fval_p,exitflag_p,output_p,residuals_p] = paretosearch(@coste_mo, 10, A, b, Aeq, beq, [], [], [], options_mo);

%% FUNCIONES 

% Coste monobjetivo
function coste_total = coste_monoobjetivo(x)

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")

M = x(1);
M_t = x(2);
AOA = x(3);
AOA_t = x(4);
Sweep = x(5);
Wingspan = x(6);
hft = x(7);
hft_t = x(8);
Fm = x(9);
E = x(10);

% Valores

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

% Aerodinámica
[CL, ~, CD, Lift, ~, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
t = 0.07*c;     % Espesor del ala

[~, ~, ~, Lift_t, Drag_t, F_beam_t, c_t, ~, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
t_t = 0.07*c_t;     % Espesor del ala

% Estructura
u = structure(E, Wingspan, c, t, F_beam);
u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);          

% Actuaciones y propulsión
[R, Aut, ~, ~, ~, ~, ~, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm*1E3, v_t, Lift_t, Drag_t, Lift, valores_T,...
                                                   valores_SFC, valores_T_af, valores_SFC_af);

% Costes
[coste_motor, ~, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, To_min, Fm*1E3, valores_E,...
                                                             precio_materiales, E, S, precio_combustible);
                    
salidas = [Aut R SFC coste_motor, coste_mat, coste_comb];
disp('    Autonomía   Alcance   SFC')
disp(['    ' num2str(salidas(1)) ' h ' num2str(salidas(2)) ' km ' num2str(salidas(3)) ' kg/kN h '])
disp(' ')
disp('    Motor   Material   Combustible')
disp(['    ' num2str(salidas(4)) ' € ' num2str(salidas(5)) ' € ' num2str(salidas(6)) ' € '])
disp(' ')

end

% Coste heurístico
function coste_total = coste_heu(x)

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")

M = x(1)/1E3;
M_t = x(2);
AOA = x(3)/1E2;
AOA_t = x(4);
Sweep = x(5)/1E2;
Wingspan = x(6)/sqrt(1E5);
hft = x(7);
hft_t = x(8);
Fm = x(9);
E = x(10);

% Valores

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

% Aerodinámica
[CL, ~, CD, Lift, ~, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
t = 0.07*c;     % Espesor del ala

[~, ~, ~, Lift_t, Drag_t, F_beam_t, c_t, ~, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
t_t = 0.07*c_t;     % Espesor del ala

% Estructura :((
u = structure(E, Wingspan, c, t, F_beam);
u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);          

% Actuaciones y propulsión
[R, Aut, ~, ~, ~, ~, ~, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm*1E3, v_t, Lift_t, Drag_t, Lift, valores_T,...
                                                   valores_SFC, valores_T_af, valores_SFC_af);

% Costes
[coste_motor, ~, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, To_min, Fm*1E3, valores_E,...
                                                             precio_materiales, E, S, precio_combustible);
                 
salidas = [Aut R SFC coste_motor, coste_mat, coste_comb];
disp('    Autonomía   Alcance   SFC')
disp(['    ' num2str(salidas(1)) ' h ' num2str(salidas(2)) ' km ' num2str(salidas(3)) ' kg/kN h '])
disp(' ')
disp('    Motor   Material   Combustible')
disp(['    ' num2str(salidas(4)) ' € ' num2str(salidas(5)) ' € ' num2str(salidas(6)) ' € '])
disp(' ')

end

% Coste multiobjetivo

function Vector_moo = coste_mo(x)

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")

M = x(1)/1E3;
M_t = x(2);
AOA = x(3)/1E2;
AOA_t = x(4);
Sweep = x(5)/1E2;
Wingspan = x(6)/sqrt(1E5);
hft = x(7);
hft_t = x(8);
Fm = x(9);
E = x(10);

% Valores

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

% Aerodinámica
[CL, ~, CD, Lift, ~, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
t = 0.07*c;     % Espesor del ala

[~, ~, ~, Lift_t, Drag_t, F_beam_t, c_t, ~, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
t_t = 0.07*c_t;     % Espesor del ala

% Estructura
u = structure(E, Wingspan, c, t, F_beam);
u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);          

% Actuaciones y propulsión
[R, Aut, ~, ~, ~, ~, ~, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm*1E3, v_t, Lift_t, Drag_t, Lift, valores_T,...
                                                   valores_SFC, valores_T_af, valores_SFC_af);

% Costes
[coste_motor, ~, coste_mat, coste_comb, coste_total_mo] = costes(valores_T, valores_T_af, precio_motores, To_min, Fm*1E3, valores_E,...
                                                             precio_materiales, E, S, precio_combustible);
% Vector MOO

Vector_moo = [coste_total_mo, -M];
                 
salidas = [Aut R SFC coste_motor, coste_mat, coste_comb];
disp('    Autonomía   Alcance   SFC')
disp(['    ' num2str(salidas(1)) ' h ' num2str(salidas(2)) ' km ' num2str(salidas(3)) ' kg/kN h '])
disp(' ')
disp('    Motor   Material   Combustible')
disp(['    ' num2str(salidas(4)) ' € ' num2str(salidas(5)) ' € ' num2str(salidas(6)) ' € '])
disp(' ')

end

