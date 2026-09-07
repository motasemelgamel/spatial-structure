% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Locate the workbook folder robustly even when Run Section causes MATLAB
% to execute a temporary copy inside an Editor_* directory.
inputFileNames = { ...
    'surviving_fraction_trajectories.xlsx'; ...
    'richness_change_trajectories.xlsx'};
candidateDirectories = {};

% The active saved Editor document is the most reliable location.
try
    activeFilename = matlab.desktop.editor.getActiveFilename;
    if ~isempty(activeFilename)
        candidateDirectories{end + 1} = fileparts(activeFilename);
    end
catch
    % MATLAB may be running without the desktop Editor.
end

% Add locations resolved from the executing file, MATLAB path, and pwd.
scriptFullPath = mfilename('fullpath');
if ~isempty(scriptFullPath)
    candidateDirectories{end + 1} = fileparts(scriptFullPath);
end
resolvedScriptPath = which('plot_fig');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

% Select the first candidate that contains both required workbooks.
scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateDirectory = candidateDirectories{directoryIndex};
    containsAllFiles = all(cellfun(@(fileName) ...
        isfile(fullfile(candidateDirectory, fileName)), inputFileNames));
    if containsAllFiles
        scriptDirectory = candidateDirectory;
        break
    end
end

if isempty(scriptDirectory)
    error(['Could not locate a folder containing both input ' ...
        'workbooks. Checked: %s'], strjoin(candidateDirectories, ', '))
end

survivalFile = fullfile(scriptDirectory, inputFileNames{1});
richnessChangeFile = fullfile(scriptDirectory, inputFileNames{2});

% Require the same theta sheets in both trajectory workbooks.
survivalSheets = string(sheetnames(survivalFile));
richnessChangeSheets = string(sheetnames(richnessChangeFile));

if ~isequal(sort(survivalSheets), sort(richnessChangeSheets))
    error('The two workbooks do not contain the same set of theta sheets.')
end

% Parse and numerically sort theta values from sheet names such as theta=0.3.
thetaValues = str2double(extractAfter(survivalSheets, 'theta='));
if any(~isfinite(thetaValues))
    error('Every worksheet name must follow the form theta=<numeric value>.')
end
[thetaValues, thetaOrder] = sort(thetaValues);
thetaSheets = survivalSheets(thetaOrder);
nTheta = numel(thetaValues);

% Reconstruct the plotted summaries directly from trajectory-level data.
meanSurvivingFractionS0 = [];
stdSurvivingFractionS0 = [];
meanRelativeRichnessChange = [];
stdRelativeRichnessChange = [];
betaValues = [];
nTrajectories = [];

for g = 1:nTheta
    sheetName = thetaSheets(g);

    [survivalBeta, survivalTrajectories] = readTrajectorySheet( ...
        survivalFile, sheetName);
    [richnessBeta, richnessTrajectories] = readTrajectorySheet( ...
        richnessChangeFile, sheetName);

    % Initialize output arrays after discovering dimensions from the files.
    if g == 1
        betaValues = survivalBeta;
        nBeta = numel(betaValues);
        nTrajectories = size(survivalTrajectories, 2);
        meanSurvivingFractionS0 = nan(nBeta, nTheta);
        stdSurvivingFractionS0 = nan(nBeta, nTheta);
        meanRelativeRichnessChange = nan(nBeta, nTheta);
        stdRelativeRichnessChange = nan(nBeta, nTheta);
    end

    assertSameBeta(betaValues, survivalBeta, ...
        sprintf('%s in surviving-fraction data', sheetName))
    assertSameBeta(betaValues, richnessBeta, ...
        sprintf('%s in richness-change data', sheetName))

    if size(survivalTrajectories, 2) ~= nTrajectories || ...
            size(richnessTrajectories, 2) ~= nTrajectories
        error('Trajectory count is inconsistent in sheet %s.', sheetName)
    end

    meanSurvivingFractionS0(:, g) = ...
        mean(survivalTrajectories, 2, 'omitnan');
    stdSurvivingFractionS0(:, g) = ...
        std(survivalTrajectories, 0, 2, 'omitnan');
    meanRelativeRichnessChange(:, g) = ...
        mean(richnessTrajectories, 2, 'omitnan');
    stdRelativeRichnessChange(:, g) = ...
        std(richnessTrajectories, 0, 2, 'omitnan');
end

fprintf(['Read %d theta sheets, %d beta values, and %d trajectories ' ...
    'per stochastic condition from %s.\n'], ...
    nTheta, numel(betaValues), nTrajectories, scriptDirectory)

% Apply the same figure formatting as optimal_beta_final.m.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
colors = lines(nTheta);
legendLabels = compose('$\\theta = %.2g$', thetaValues);
richnessMarkerSize = 15;
richnessLineWidth = 3;
richnessFontSize = 75;

leftPanelMinimum = min( ...
    meanSurvivingFractionS0 - stdSurvivingFractionS0, ...
    [], 'all', 'omitnan');
rightPanelMinimum = min( ...
    meanRelativeRichnessChange - stdRelativeRichnessChange, ...
    [], 'all', 'omitnan');
if ~isfinite(leftPanelMinimum) || leftPanelMinimum >= 1
    leftPanelMinimum = 0.9;
end
if ~isfinite(rightPanelMinimum) || rightPanelMinimum >= 0
    rightPanelMinimum = -0.1;
end

% Figure 1: stochastic richness metrics shown as two square panels.
figure(1)
clf
set(gcf, 'Color', 'w', 'Position', [50 80 1900 900])
tiledlayout(1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: fraction surviving relative to the original species count.
nexttile
hold on
richnessHandles = gobjects(nTheta, 1);
for g = 1:nTheta
    richnessHandles(g) = errorbar( ...
        betaValues, meanSurvivingFractionS0(:, g), ...
        stdSurvivingFractionS0(:, g), ...
        'Color', colors(g, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', richnessMarkerSize, ...
        'MarkerFaceColor', colors(g, :), ...
        'MarkerEdgeColor', colors(g, :), ...
        'LineWidth', richnessLineWidth, 'CapSize', 10);
end
formatRichnessPanel( ...
    'Fraction of surviving species', richnessFontSize, betaValues, ...
    [leftPanelMinimum 1])
legend(richnessHandles, legendLabels, ...
    'Location', 'best', 'Orientation', 'vertical', 'FontSize', 42)

% Panel 2: fractional richness change relative to the burn-in state.
nexttile
hold on
for g = 1:nTheta
    errorbar(betaValues, meanRelativeRichnessChange(:, g), ...
        stdRelativeRichnessChange(:, g), ...
        'Color', colors(g, :), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', richnessMarkerSize, ...
        'MarkerFaceColor', colors(g, :), ...
        'MarkerEdgeColor', colors(g, :), ...
        'LineWidth', richnessLineWidth, 'CapSize', 10);
end
formatRichnessPanel( ...
    'Fractional change of species from steady state', ...
    richnessFontSize, betaValues, [rightPanelMinimum 0])

function [betaValues, trajectoryValues] = ...
    readTrajectorySheet(workbookPath, sheetName)

% Read one beta column followed by all available trajectory columns.
dataTable = readtable(workbookPath, 'Sheet', char(sheetName), ...
    'VariableNamingRule', 'preserve');
variableNames = string(dataTable.Properties.VariableNames);
if width(dataTable) < 2 || ~strcmpi(variableNames(1), 'beta')
    error('Sheet %s in %s must start with beta and trajectory columns.', ...
        sheetName, workbookPath)
end
if any(~startsWith(variableNames(2:end), 'trajectory_'))
    error('Sheet %s in %s contains unexpected trajectory headers.', ...
        sheetName, workbookPath)
end

betaValues = dataTable{:, 1};
trajectoryValues = dataTable{:, 2:end};
if ~isnumeric(betaValues) || ~isnumeric(trajectoryValues)
    error('Sheet %s in %s must contain numeric data.', ...
        sheetName, workbookPath)
end
betaValues = betaValues(:);
end

function assertSameBeta(referenceBeta, candidateBeta, contextText)

% Require identical beta grids up to floating-point roundoff.
candidateBeta = candidateBeta(:);
if numel(candidateBeta) ~= numel(referenceBeta)
    error('Beta-grid length mismatch for %s.', contextText)
end
comparisonScale = max(1, max(abs(referenceBeta), [], 'omitnan'));
tolerance = 100 .* eps(comparisonScale);
if any(abs(candidateBeta - referenceBeta) > tolerance | ...
        xor(isnan(candidateBeta), isnan(referenceBeta)))
    error('Beta values do not match for %s.', contextText)
end
end

function formatRichnessPanel(yLabelText, fontSize, betaValues, yLimits)

% Apply the shared scientific formatting used by both richness panels.
set(gca, 'XScale', 'log', 'FontSize', fontSize, 'LineWidth', 1.5)
xlabel('$\beta$')
ylabel(yLabelText)
xlim([min(betaValues) max(betaValues)])
ylim(yLimits)
xticks([1e-2 1 1e2])
pbaspect([1 1 1])
box on
end
