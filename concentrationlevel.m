% File: nanoparticle_dynamics.m

% Clear workspace and figures
clear;
close all;

% Define parameters (Change these values based on your data)
Kmaxi = 1; % Maximum uptake rate constant
t50i = 1; % Time to reach half-maximum uptake rate
ni = 1; % Hill coefficient
Qi = 1; % Regional blood flow
Ca = 1; % Concentration of NMs in arterial plasma
PAi = 1; % Product of permeability coefficient
Pi = 1; % Tissue-plasma distribution coefficient
Krei = 1; % Release rate constant

% Time vector for simulation
tspan = [0 24]; % 24 hours simulation

% Initial conditions [ABloodi, ATi, APCi_or_ATCi]
initial_conditions = [0, 0, 0]; % Change as necessary

% Solve differential equations
[time, results] = ode45(@(t, y) nanoparticle_model(t, y, Kmaxi, t50i, ni, Qi, Ca, PAi, Pi, Krei), tspan, initial_conditions);

% Plot results
figure;
plot(time, results(:, 1), 'b-', 'LineWidth', 2);
hold on;
plot(time, results(:, 2), 'r--', 'LineWidth', 2);
plot(time, results(:, 3), 'g-.', 'LineWidth', 2);
hold off;

xlabel('Time (h)');
ylabel('Concentration (units)');
title('Nanoparticle Dynamics in Tumor-Bearing Mice');
legend('Capillary Concentration', 'Interstitium Concentration', 'Endocytic/Phagocytic Cell Concentration', 'Location', 'Best');
grid on;

% Save the figure
saveas(gcf, 'nanoparticle_dynamics.png');

% Differential equations model
function dydt = nanoparticle_model(t, y, Kmaxi, t50i, ni, Qi, Ca, PAi, Pi, Krei)
    % Unpack variables
    ABloodi = y(1);
    ATi = y(2);
    APCi_or_ATCi = y(3);

    % Calculate intermediate variables
    Kupi = Kmaxi .* (t.^ni) ./ (t50i^ni + t.^ni);
    CVi = ABloodi; % Assuming CVi approximates to ABloodi for simplification
    CTi = ATi; % Assuming CTi approximates to ATi for simplification

    % Equations
    dABloodi_dt = Qi * (Ca - CVi) + PAi * (CTi - CVi / Pi);
    dATi_dt = PAi * CVi - (PAi / Pi) * CTi - Kupi * ATi + Krei * APCi_or_ATCi;
    dAPCi_or_ATCi_dt = Kupi * ATi - Krei * APCi_or_ATCi;

    dydt = [dABloodi_dt; dATi_dt; dAPCi_or_ATCi_dt];
end