%% D-CELL DESIIGN =========================================================
clear
clc
close all

b = 65;
fuselage_diam = 6.49;   % fuselage diameter (m)
fuselage_rad = fuselage_diam/2; % fuselage radius (m)
y = 0:0.1:65/2;
y_wing = y(33:end);
c_root = 13.6;
c_tip = 2.16;

c = linspace(c_root,c_tip,length(y_wing));      % Wing chord (m)
b_ref = b/2 - fuselage_rad;   % Wingbox span
r_DCell = 0.004 .* c;        % D-Cell radius
b_curve = D_Cell_Arc(c);    % D-Cell arc length
U_cruise = 0.83 * 340;
E = 72e9; 

rho_Al = 2700;

n_psrib = 10:2:50;               % Number of Pseudo Ribs
t_DCell = 1.6:0.05:2.2;              % D-Cell thickness (mm)
t_rib = 0.5; %mm


% Itterate No Pseudo Ribs and Thickness
for n = 1:length(n_psrib)

    a = b_ref/(n_psrib(n) + 1);    % D-Cell span length

    for t = 1:length(t_DCell)
        
        % Determine Ks
        Ks = curved_Ks(a,b_curve,r_DCell,t_DCell(t));
        %Ks = ones(size(b_curve))*20;
        
        Buckling_stress(n,t) = max(Ks .* E .* ((t_DCell(t)./1000)./b_curve).^2 ./ (10^6));


        % Rib Mass
        c_rib = linspace(c_root,c_tip,n_psrib(n));
        A = 0.0166 .* c_rib.^2;
        m_rib(n) = sum((t_rib/1000 .* A) .* rho_Al);

        % Skin Mass
        m_skin(t) = sum((t_DCell(t)/1000) .* b_curve .* 0.1 .* rho_Al) ;

        % Total Mass
        m_total(n,t) = m_rib(n) + m_skin(t);
       


    end
end

Buckling_stress
m_total

% Plot Stress
figure
hold on
surf(t_DCell,n_psrib,Buckling_stress)
%shading interp  % Smooth shading
Z = ones(size(Buckling_stress)) * 4.7; % Create a flat plane at Z = 4.7
surf(t_DCell,n_psrib, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'k')

xlabel('Number of Pseudo Ribs')
ylabel('D-Cell Thickness (mm)')
zlabel('Buckling Stress (MPa)')
title('Buckling Stress Distribution')
view(3)
colorbar


% % Plot Mass Opt
% figure
% hold on
% surf(t_DCell,n_psrib,m_total)
% shading interp  % Smooth shading
% xlabel('Number of Pseudo Ribs')
% ylabel('D-Cell Thickness (mm)')
% zlabel('Mass (Kg)')
% title('Mass Optimisation')
% colorbar


% Compute efficiency metric: Higher value means better trade-off
efficiency = Buckling_stress ./ m_total;

% Find the best combination of ribs and skin thickness
[max_eff, idx] = max(efficiency(:));
[row, col] = ind2sub(size(efficiency), idx);
best_n_psrib = n_psrib(row);
best_t_DCell = t_DCell(col);

fprintf('Best Combination:\n');
fprintf('Number of Pseudo Ribs: %d\n', best_n_psrib);
fprintf('D-Cell Thickness: %.2f mm\n', best_t_DCell);
fprintf('Highest Efficiency: %.4f MPa/kg\n', max_eff);

% Plot Efficiency Surface
figure
surf(t_DCell, n_psrib, efficiency)
%shading interp
xlabel('D-Cell Thickness (mm)')
ylabel('Number of Pseudo Ribs')
zlabel('Efficiency (MPa/kg)')
title('Efficiency: Buckling Stress per Unit Mass')
colorbar
view(3)
hold on

% % Highlight Best Point
% plot3(best_t_DCell, best_n_psrib, max_eff, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
% text(best_t_DCell, best_n_psrib, max_eff, sprintf('Best: %d ribs, %.2f mm', best_n_psrib, best_t_DCell), ...
%     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r')


% figure
% plot(n_psrib,m_rib)
% figure
% plot(t_DCell,m_skin)










%% === FUNCTIONS === %%

function arc = D_Cell_Arc(c)

    % Load the airfoil data
    filename = 'sc20712.dat.txt';
    data = load(filename);
    
    % Extract x and y coordinates
    x = data(:,1);
    y = data(:,2);

    % Define the indices for the leading edge and wing box boundary
    start_idx = 23;
    end_idx = 184;
    
    x2 = [x(start_idx:-1:1); x(end:-1:end_idx)];
    y2 = [y(start_idx:-1:1); y(end:-1:end_idx)];
    
    % Compute the cumulative distance along the airfoil
    arc_length = 0;
    for i = 1:length(x2)-1
        dx = x2(i+1) - x2(i);
        dy = y2(i+1) - y2(i);
        arc_length = arc_length + sqrt(dx^2 + dy^2);
    end

    arc = arc_length .* c;

end

function Ks = curved_Ks(a,b_curve,r_DCell,t_DCell)
    
    for i = 1:length(b_curve)
        %a_b = a./b_curve
        if a > b_curve(i)
        %if a_b >= 1
            a_b(i) = a./b_curve(i);
            b_rt(i) = b_curve(i) ./ sqrt(r_DCell(i).*t_DCell);
            %a_rt(i) = 0;

            if 14 < b_rt(i) && b_rt(i) < 15
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 44.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 41;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 37;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 33;
                elseif 3 < a_b(i)
                    Ks(i) = 27.5;
                end

            elseif 13 < b_rt(i) && b_rt(i) < 14
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 41;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 37;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 33.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 30;
                elseif 3 < a_b(i)
                    Ks(i) = 25.5;
                end

            elseif 12 < b_rt(i) && b_rt(i) < 13
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 36.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 33.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 30;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 27.5;
                elseif 3 < a_b(i)
                    Ks(i) = 23;
                end

            elseif 11 < b_rt(i) && b_rt(i) < 12
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 33;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 30.1;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 27;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 24.5;
                elseif 3 < a_b(i)
                    Ks(i) = 21;
                end

            elseif 10 < b_rt(i) && b_rt(i) < 11
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 29.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 27;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 24.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 22;
                elseif 3 < a_b(i)
                    Ks(i) = 19;
                end
    
            elseif 9 < b_rt(i) && b_rt(i) < 10
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 26.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 24;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 21.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 20;
                elseif 3 < a_b(i)
                    Ks(i) = 17;
                end
            
            elseif 8 < b_rt(i) && b_rt(i) < 9
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 23.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 21.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 18.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 17;
                elseif 3 < a_b(i)
                    Ks(i) = 16;
                end
    
            elseif 7 < b_rt(i) && b_rt(i) < 8
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 20.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 18;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 16.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 15;
                elseif 3 < a_b(i)
                    Ks(i) = 13.5;
                end
    
            elseif 6 < b_rt(i) && b_rt(i) < 7
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 17;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 15.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 13.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 12.5;
                elseif 3 < a_b(i)
                    Ks(i) = 11;
                end
    
            elseif 5 < b_rt(i) && b_rt(i) < 6
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 14.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 13;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 11.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 10.5;
                elseif 3 < a_b(i)
                    Ks(i) = 9.5;
                end
    
            elseif 3 < b_rt(i) && b_rt(i) < 5
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 11;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 10;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 8.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 8;
                elseif 3 < a_b(i)
                    Ks(i) = 7.5;
                end
    
            elseif b_rt(i) < 3 && b_rt(i) ~= 0
    
                if  a_b(i) < 1.25 && a_b(i) ~= 0
                    Ks(i) = 8.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 7.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 6.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 6;
                elseif 3 < a_b(i)
                    Ks(i) = 5.5;
                end

            elseif b_rt >15
                disp('b/rt out of bounds (high)')
            else
                %disp('b/rt zero')
            end
        

        % ================================================================

        elseif b_curve(i) > a
            b_a(i) = b_curve(i)./a;
            a_rt(i) = a ./ sqrt(r_DCell(i).*t_DCell);

            if 17 < a_rt(i) && a_rt(i) < 18
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 58.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 55.5;
                else
                    Ks(i) = 54;
                end

            elseif 16 < a_rt(i) && a_rt(i) < 17
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 54.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 51;
                else
                    Ks(i) = 49.5;
                end

            elseif 15 < a_rt(i) && a_rt(i) < 16
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 50;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 46.5;
                else
                    Ks(i) = 45;
                end

            elseif 14 < a_rt(i) && a_rt(i) < 15
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 45.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 42;
                else
                    Ks(i) = 40.5;
                end

            elseif 13 < a_rt(i) && a_rt(i) < 14
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 41;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 38;
                else
                    Ks(i) = 36.5;
                end

            elseif 12 < a_rt(i) && a_rt(i) < 13
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 37;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 34;
                else
                    Ks(i) = 32.5;
                end

            elseif 11 < a_rt(i) && a_rt(i) < 12
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 33.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 30.25;
                else
                    Ks(i) = 29;
                end

            elseif 10 < a_rt(i) && a_rt(i) < 11
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 30;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 27;
                else
                    Ks(i) = 25.5;
                end
    
            elseif 9 < a_rt(i) && a_rt(i) < 10
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 26.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 23.5;
                else
                    Ks(i) = 22;
                end

            elseif 8 < a_rt(i) && a_rt(i) < 9
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 23;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 20;
                else
                    Ks(i) = 18.5;
                end
    
            elseif 7 < a_rt(i) && a_rt(i) < 8
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 20;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 17.5;
                else
                    Ks(i) = 16;
                end

            elseif 6 < a_rt(i) && a_rt(i) < 7
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 17;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 14.5;
                else
                    Ks(i) = 13;
                end

            elseif 5 < a_rt(i) && a_rt(i) < 6
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 14.5;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 12;
                else
                    Ks(i) = 10.75;
                end

            elseif 4 < a_rt(i) && a_rt(i) < 5
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 12;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 9;
                else
                    Ks(i) = 8;
                end

            elseif 3 < a_rt(i) && a_rt(i) < 4
    
                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 10;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 8;
                else
                    Ks(i) = 7;
                end
    
            elseif 2 < a_rt(i) && a_rt(i) < 3

                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 9;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 7.2;
                else
                    Ks(i) = 6.3;
                end

            elseif 0 < a_rt(i) && a_rt(i) < 2

                if  b_a(i) < 1.5 && b_a(i) ~= 0
                    Ks(i) = 7.8;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 6.5;
                else
                    Ks(i) = 5.5;
                end


            elseif a_rt ==0
                disp('a/rt zero')

            else
                %disp('a/rt out of bounds (high)')
            end
        end
    end
end