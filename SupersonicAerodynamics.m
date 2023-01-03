function [C_L, C_M, C_D, Lift, Drag, F_beam, c, WS, v] = SupersonicAerodynamics(M, AOA, Sweep, b, hft)

    DOF = 202;
    [~, a, ~, rho] = atmosisa(hft/3.28084);
    beta = sqrt(M^2-1);
    v = M*a;
    alfa_rad = AOA * (pi / 180);
    sweep = acotd(cotd(Sweep)/beta);
    L = 1:80;
    N = -50:50;

    for i = 1:max(L) 
        for j = 0:max(N)
            x_le(j+1) = j / cotd(sweep);
            x_te(j+1) = max(N) / cotd(sweep);


            if L(i) - x_le(j+1) <= 0 
                A(i, j+1) = 0;
            elseif L(i) - x_le(j+1) >= 1 
                A(i, j+1) = 1;
            else
                A(i, j+1) = L(i) - x_le(j+1); 
            end

            if L(i) - x_te(j+1) >= 1 
                B(i, j+1) = 0;
            elseif L(i) - x_te(j+1) <= 0 
                B(i, j+1) = 1;
            else
                B(i, j+1) = 1 - (L(i) - x_te(j+1));
            end

            C = ones(max(L), max(N)+1);
            C(max(L), max(N)+1) = 0.5;

        end
    end

    A = [fliplr(A), A(:, 2:end)];
    B = [fliplr(B), B(:, 2:end)];
    C = [fliplr(C), C(:, 2:end)];
    x_le = [flip(x_le), x_le(2:end)];
    x_te = [flip(x_te), x_te(2:end)];
    Cp = zeros(length(L), length(N));
    Cp_b = zeros(length(L), length(N)); 
    
    for i = 1:max(L)
        for j = max(N):(length(N) - 1)
            if and((L(i) - x_le(j+1)) >= 0, round(L(i)-x_te(j+1)) <= 0) 
                R_a = zeros(length(L), length(N));
                for m = 1:(i - 1)
                    for n = 1:length(N)
                        if and(abs(atand((N(n) - N(j+1))/(L(m) - L(i)))) <= 45, L(m)-x_le(n) >= 0)
                            R_a(m, n) = ((L(i) - L(m) + 0.5)^2 - ...
                                (N(j+1) - N(n) - 0.5)^2)^0.5 ...
                                / ((L(i) - L(m) + 0.5) * (N(j+1) - N(n) ...
                                -0.5)) - ((L(i) - L(m) + 0.5)^2 - ...
                                (N(j+1) - N(n) + 0.5)^2)^0.5 / ((L(i) - ...
                                L(m) + 0.5) * (N(j+1) - N(n) + 0.5));

                        else
                            R_a(m, n) = 0; 
                        end
                    end
                end
                Cp(i, j+1) = 4 / beta * alfa_rad + 1 / pi * sum(sum(R_a.*A.*B.*C.*Cp));
                Cp = [fliplr(Cp(:, (max(N) + 2):end)), Cp(:, (max(N) + 1):end)];

            else
                Cp(i, j+1) = 0;
            end
        end
        
        for j = max(N):(length(N) - 1)
            if i ~= length(L)
                R_b = zeros(length(L), length(N)); 
                for m = 1:i
                    for n = 1:length(N)
                        if and(abs(atand((N(n) - N(j+1))/(L(m) - L(i+1)))) <= 45, ...
                                L(m)-x_le(n) > 0)

                            R_b(m, n) = ((L(i+1) - L(m) + 0.5)^2 - (N(j+1) - N(n) ...
                                -0.5)^2)^0.5 / ((L(i+1) - L(m) + 0.5) * ...
                                (N(j+1) - N(n) - 0.5)) - ((L(i+1) - L(m) + 0.5)^2 - (N(j+1) - N(n) + 0.5)^2)^0.5 / ...
                                ((L(i+1) - L(m) + 0.5) * (N(j+1) - N(n) + 0.5));
                        else
                            R_b(m, n) = 0;

                        end
                    end
                end
            else
                R_b(i, j+1) = 0;
            end
            if or(L(i)-x_le(j+1) < 0, round(L(i)-x_te(j+1)) > 0) 
                Cp_b(i, j+1) = 0;
            else
                Cp_b(i, j+1) = 4 / beta * alfa_rad + 1 / pi * sum(sum(R_b.*A.*B.*C.*Cp));
                Cp_b = [fliplr(Cp_b(:, (max(N) + 2):end)), Cp_b(:, (max(N) + 1):end)];
            end
        end
        
        for j = max(N):(length(N) - 1)

            if L(i) - x_le(j+1) <= 1
                Cp(i, j+1) = (1 / 2) * (1 + A(i, j+1) / (1 + A(i, j+1))) * Cp(i, j+1) + (1 / 2) * ...
                    (A(i, j+1) / (1 + A(i, j+1))) * Cp_b(i, j+1); 
            else
                Cp(i, j+1) = (3 / 4) * Cp(i, j+1) + (1 / 4) * Cp_b(i, j+1); 
            end
            Cp = [fliplr(Cp(:, (max(N) + 2):end)), Cp(:, (max(N) + 1):end)]; 
        end
    end
    %Wing Area
    S_ = (2 / beta) * sum(sum(A.*B.*C));
    %Wing span
    b_ = length(N);
    %Wing Mean Chord
    c_ = S_ / b_;
    %Wing Slope
    slope = zeros(max(L), length(N));
    slope(:, :) = -alfa_rad;
    slope_L = [slope(2:end, :); zeros(1, length(N))];
    %Cp(L+1,N)
    Cp_L = [Cp(2:end, :); zeros(1, length(N))];

    %AERODYNAMIC COEFFICIENTS FOR THE WING
    %Lift Coefficient
    C_L = (2 / (beta * S_)) * sum(sum(((3 / 4) * Cp + (1 / 4) * Cp_L).*A.*B.*C));
    %Pitching-moment Coefficient
    C_M = (2 / (beta * S_ * c_)) * sum(sum(transpose(L).*((3 / 4) * Cp + (1 / 4) * Cp_L).*A.*B.*C));
    %Drag Coefficient
    C_D = 0.018 - (2 / (beta * S_)) * sum(sum(((3 / 4) * Cp + (1 / 4) * Cp_L).*((3 / 4) * slope + (1 / 4) * slope_L).*A.*B.*C));

    %REDIMESIONALIZATION
    % Wing surface
    WS = b^2*tand(Sweep)/4;
    c = WS/b;
    % Lift
    Lift = 1/2*rho*v^2*WS*C_L;
    % Drag
    Drag = 1/2*rho*v^2*WS*C_D;

    %EQUIVALENT FORCE
    F_beam = zeros(DOF,1);
    F_beam(round(DOF/3),:) = Lift;

end
