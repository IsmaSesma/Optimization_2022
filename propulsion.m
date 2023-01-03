%% Función propulsión

function [SFC] = propulsion(valores_T, valores_T_af, valores_SFC, valores_SFC_af, T)
    
    % valores_T - vector con los datos de empuje de los motores [N]
    % valores_T_af - vector con los datos de empuje CON POSTCOMBUSTOR de los motores [N]
    % valores_SFC - vector con los datos de consumo específico de los motores [kg/kNh]
    % valores_SFC_af - vector con los datos de consumo específico CON POSTCOMBUSTOR de los motores [kg/kNh]
    % T - empuje [N]

    % si el valor de T es mayor que los de empuje de los motores sin PC se va a los empujes con PC

    if T < max(valores_T)

        x = linspace(0, max( valores_T ));
    
        coefs = polyfit( valores_T, valores_SFC, 2 ); % coefs polinomio para aproximar
        
        funcion_SFC = coefs(1)*x.^2 + coefs(2)*x + coefs(3);
        SFC = coefs(1)*T^2 + coefs(2)*T + coefs(3);

        disp('No hace falta postcombustor')

%         figure
%         plot(valores_T, valores_SFC, '*')
%         hold on
%         plot(x, funcion_SFC, 'LineWidth',1)
%         title('Sin postcombustor')
%         xlabel('T [N]')
%         ylabel('SFC [kg/kNh]')
%         hold off

    else

        x = linspace(0, max( valores_T_af ));
    
        coefs = polyfit( valores_T_af, valores_SFC_af, 2 ); % coefs polinomio para aproximar
        
        funcion_SFC = coefs(1)*x.^2 + coefs(2)*x + coefs(3);
        SFC = coefs(1)*T^2 + coefs(2)*T + coefs(3);

        disp('Es necesario postcombustor')

%         figure
%         plot(valores_T_af, valores_SFC_af, '*')
%         hold on
%         plot(x, funcion_SFC, 'LineWidth',1)
%         title('Con postcombustor')
%         xlabel('T [N]')
%         ylabel('SFC [kg/kNh]')
%         hold off

    end

end