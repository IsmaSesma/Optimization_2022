% Función de costes

function [coste_motor, precio_mat, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, T, Fm, ...
                                                                                valores_E, precio_materiales, E, S, precio_combustible)

    % valores_T - vector con los datos de empuje de los motores [N]
    % precio_motores - vector precio de los motores [€]
    % T - empuje [N]
    % valores_E - vector con los datos de modulo elástico de los materiales [Pa]
    % precio_materiales - vector con los precios de los materiales [€/m^3]
    % E - modulo elástico [Pa]
    % t - espesor MISMAS UNIDADES QUE LA SUPERFICIE [m]
    % S- superficie alar [m^2]
    % precio_combustible - precio del combustible JET-B [€/L]
    % SFC - consumo específico [kg/kNh]


    % Coste motor

    if T < max(valores_T)

        x_T = linspace(0, max( valores_T ));
    
        coefs = polyfit( valores_T, precio_motores, 2 ); % coefs polinomio para aproximar
        pl = polyval( coefs, x_T ); % valores del polinomio en cada x
        
        funcion_T =  coefs(1)*x_T.^2 + coefs(2)*x_T + coefs(3);
        coste_motor = coefs(1)*T^2 + coefs(2)*T + coefs(3);
     
%         figure
%         plot(valores_T, precio_motores, '*')
%         hold on
%         plot(x_T, funcion_T, 'LineWidth',1)
%         title('Sin postcombustor')
%         xlabel('T [N]')
%         ylabel('Precio')
%         hold off

    else

        x_T = linspace(0, max(valores_T_af));
    
        coefs = polyfit(valores_T_af, precio_motores, 1); % coefs polinomio para aproximar
        pl = polyval( coefs, x_T ); % valores del polinomio en cada x
        
        funcion_T =  coefs(1)*x_T + coefs(2);
        coste_motor = coefs(1)*T + coefs(2);
    
%         figure
%         plot(valores_T, precio_motores, '*')
%         hold on
%         plot(x_T, funcion_T, 'LineWidth',1)
%         title('Con postcombustor')
%         xlabel('T [N]')
%         ylabel('Precio')
%         hold off

    end

    % Coste material

    % 1- Ti6Al4V
    % 2- 7068 Al
    % 3- 6061 Al
    % 4- invar
    % 5- duraluminio

    E = E*1E9;
    x_mat = linspace(0, max(valores_E));

    coefs = polyfit(valores_E, precio_materiales, 2); % coefs polinomio para aproximar
    pl = polyval( coefs, x_mat ); % valores del polinomio en cada x
    
    funcion_mat = coefs(1)*x_mat.^2 + coefs(2)*x_mat + coefs(3);
    precio_mat = coefs(1)*E^2 + coefs(2)*E + coefs(3); % €/m^3

    Vol = 2*S*0.02; % m^3
    coste_mat = Vol*precio_mat; % €

%     figure
%     plot(valores_E, precio_materiales, '*')
%     hold on
%     plot(x_mat, funcion_mat, 'LineWidth',1)
%     xlabel('E [Pa]')
%     ylabel('Precio')
%     hold off

    % Coste combustible
    
    rho_JetB = 800;     % kg/m^3
    Vol = Fm/rho_JetB*1E3;
    coste_comb = Vol*precio_combustible; 

    % Coste total

    coste_total = coste_motor + coste_mat + coste_comb; % €

end