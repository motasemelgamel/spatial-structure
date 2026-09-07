% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all

% Define the workbook locations and parameter values represented by their tabs.
scriptFile = matlab.desktop.editor.getActiveFilename;
dataDirectory = fileparts(scriptFile);
betaSweepWorkbook = fullfile(dataDirectory, "beta_sweeps_at_fixed_theta.xlsx");
thetaSweepWorkbook = fullfile(dataDirectory, "theta_sweeps_at_fixed_beta.xlsx");
abundanceWorkbook = fullfile(dataDirectory, "species_abundances_by_theta.xlsx");

betaValues = [0.001 0.01 0.1 1 10 100 1000];
thetaValues = [0 0.3 0.5 0.8 1];
abundanceThetaValues = [0 0.5 1];
numberOfSpecies = 50;

plotColors = [0.27 0.15 0.28; ...
              0.44 0.53 0.46; ...
              0.88 0.00 0.00; ...
              0.64 0.96 0.24; ...
              0.68 0.29 0.67; ...
              0.70 0.07 0.25; ...
              0.22 0.67 0.84];

% Plot the surviving fraction across theta for every fixed-beta sheet.
figure('Name', 'Survival versus theta', 'NumberTitle', 'off')
hold on
for betaIndex = 1:numel(betaValues)
    sheetName = "beta=" + betaValues(betaIndex);
    data = readtable(thetaSweepWorkbook, 'Sheet', sheetName);
    plot(data.theta, data.surviving_fraction, ...
        'Color', plotColors(betaIndex, :), ...
        'Marker', '.', 'LineWidth', 3, ...
        'LineStyle', 'none', 'MarkerSize', 45)
end
legend("$\beta=" + betaValues + "$", 'Interpreter', 'latex')
xlabel('$\theta$', 'Interpreter', 'latex')
ylabel('Fraction of surviving species', 'Interpreter', 'latex')
ylim([0 1])
set(gca, 'FontSize', 65)
pbaspect([1 1 1])
box on
hold off

% Plot the surviving fraction across beta for every fixed-theta sheet.
figure('Name', 'Survival versus beta', 'NumberTitle', 'off')
hold on
for thetaIndex = 1:numel(thetaValues)
    sheetName = "theta=" + thetaValues(thetaIndex);
    data = readtable(betaSweepWorkbook, 'Sheet', sheetName);
    plot(data.beta, data.surviving_fraction, ...
        'Color', plotColors(thetaIndex, :), ...
        'Marker', '.', 'LineWidth', 3, ...
        'LineStyle', 'none', 'MarkerSize', 45)
end
legend("$\theta=" + thetaValues + "$", 'Interpreter', 'latex')
xlabel('$\beta$', 'Interpreter', 'latex')
ylabel('Fraction of surviving species', 'Interpreter', 'latex')
ylim([0 1])
xticks([1e-2 1 1e2])
set(gca, 'FontSize', 70, 'XScale', 'log')
pbaspect([1 1 1])
box on
hold off

% Reuse seeded species colors so corresponding species match across theta figures.
rng(1)
speciesColors = rand(numberOfSpecies, 3);

% Plot each abundance sheet in its own figure while overlaying all species.
for thetaIndex = 1:numel(abundanceThetaValues)
    thetaValue = abundanceThetaValues(thetaIndex);
    sheetName = "theta=" + thetaValue;
    data = readtable(abundanceWorkbook, 'Sheet', sheetName);
    time = data.time;
    abundances = data{:, 2:end};
    maximumAbundance = max(abundances, [], 'all');
    abundances(abundances < 1e-3) = 0;

    figure('Name', "Abundances, theta=" + thetaValue, 'NumberTitle', 'off')
    hold on
    for speciesIndex = 1:size(abundances, 2)
        plot(time, abundances(:, speciesIndex), ...
            'Color', speciesColors(speciesIndex, :), ...
            'LineWidth', 3, 'MarkerSize', 45)
    end
    yticks([1e-2 1e0 1e2])
    ylim([1e-3 maximumAbundance])
    xlim([0 1e3])
    xlabel('Time', 'Interpreter', 'latex')
    ylabel('Population size', 'Interpreter', 'latex')
    set(gca, 'FontSize', 65, 'YScale', 'log')
    title("$\theta=" + thetaValue + "$", 'Interpreter', 'latex')
    pbaspect([1 1 1])
    box on
    hold off
end
