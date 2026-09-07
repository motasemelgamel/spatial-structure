% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all

% Define the Figure 4 workbooks and their scenario tabs.
scriptFile = matlab.desktop.editor.getActiveFilename;
dataDirectory = fileparts(scriptFile);
thetaScenarioWorkbook = fullfile(dataDirectory, "survival_vs_theta_by_nutrient_scenario_different_preferences.xlsx");
poorNutrientWorkbook = fullfile(dataDirectory, "survival_vs_number_of_poor_nutrients_different_preferences.xlsx");

scenarioSheetNames = ["alpha_only", ...
                      "delta1_only", ...
                      "alpha_delta1", ...
                      "alpha_delta1_delta2", ...
                      "alpha_delta1_delta2_delta3"];
scenarioTitles = ["$\alpha$", ...
                  "$\delta$", ...
                  "$\alpha+\delta_1$", ...
                  "$\alpha+\delta_1+\delta_2$", ...
                  "$\alpha+\delta_1+\delta_2+\delta_3$"];
highlightColors = [0.27 0.15 0.28; ...
                   0.75 0.00 0.25; ...
                   0.44 0.53 0.46; ...
                   0.88 0.52 0.94; ...
                   0.64 0.96 0.24];
previousCurveColor = [0.72 0.72 0.72];

% Read each scenario once so previous curves can be reused as gray dots.
scenarioData = cell(size(scenarioSheetNames));
for scenarioIndex = 1:numel(scenarioSheetNames)
    scenarioData{scenarioIndex} = readtable(thetaScenarioWorkbook, ...
        'Sheet', scenarioSheetNames(scenarioIndex));
end

% The alpha-containing scenarios correspond to zero through three poor nutrients.
alphaCumulativeIndices = [1 3 4 5];

% Plot each scenario separately and gray out every lower-poor-nutrient curve.
for scenarioIndex = 1:numel(scenarioSheetNames)
    figure('Name', scenarioSheetNames(scenarioIndex), ...
        'NumberTitle', 'off', 'Color', 'w', ...
        'Position', [100 + 40*scenarioIndex, 100 + 40*scenarioIndex, 900, 900])
    hold on

    cumulativePosition = find(alphaCumulativeIndices == scenarioIndex, 1);
    if ~isempty(cumulativePosition) && cumulativePosition > 1
        previousIndices = alphaCumulativeIndices(1:cumulativePosition - 1);
        for previousIndex = previousIndices
            previousData = scenarioData{previousIndex};
            plot(previousData.theta, previousData.surviving_fraction, ...
                'Color', previousCurveColor, ...
                'Marker', '.', 'LineWidth', 3, ...
                'LineStyle', 'none', 'MarkerSize', 45)
        end
    end

    currentData = scenarioData{scenarioIndex};
    plot(currentData.theta, currentData.surviving_fraction, ...
        'Color', highlightColors(scenarioIndex, :), ...
        'Marker', '.', 'LineWidth', 3, ...
        'LineStyle', 'none', 'MarkerSize', 45)
    xlabel('$\theta$', 'Interpreter', 'latex')
    ylabel({'Fraction of', 'surviving species'}, 'Interpreter', 'latex')
    xlim([0 1])
    ylim([0 1])
    xticks([0 0.5 1])
    yticks([0 0.5 1])
    set(gca, 'FontSize', 85, 'TickLabelInterpreter', 'latex')
    title(scenarioTitles(scenarioIndex), 'Interpreter', 'latex')
    pbaspect([1 1 1])
    box on
    hold off
end

% Compare poor-nutrient counts for the three fixed-theta sheets.
thetaValues = [0 0.5 1];
thetaColors = [0 0 0; 0 0 1; 1 0 1];

figure('Name', 'Survival versus poor nutrient count', ...
    'NumberTitle', 'off', 'Color', 'w', ...
    'Position', [340, 340, 900, 900])
hold on
for thetaIndex = 1:numel(thetaValues)
    sheetName = "theta=" + thetaValues(thetaIndex);
    data = readtable(poorNutrientWorkbook, 'Sheet', sheetName);
    plot(data.number_of_poor_nutrients, data.surviving_fraction, ...
        'Color', thetaColors(thetaIndex, :), ...
        'Marker', 'o', 'MarkerFaceColor', thetaColors(thetaIndex, :), ...
        'LineWidth', 3, 'MarkerSize', 20)
end
xlabel('No. of poor nutrients', 'Interpreter', 'latex')
ylabel({'Fraction of', 'surviving species'}, 'Interpreter', 'latex')
xlim([0 5])
ylim([0 1])
xticks([0 5])
yticks([0 0.5 1])
set(gca, 'FontSize', 90, 'TickLabelInterpreter', 'latex')
legend("$\theta=" + thetaValues + "$", ...
    'Interpreter', 'latex', 'Location', 'southeast', ...
    'FontSize', 45, 'Box', 'on')
pbaspect([1 1 1])
box on
hold off
