% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Reproduce the species-versus-eta figure using only the saved Excel data.
inputFileName = 'species_vs_eta.xlsx';

% Locate the workbook beside this script, including when MATLAB runs a
% temporary Editor_* copy from another working directory.
candidateDirectories = {};
try
    activeFilename = matlab.desktop.editor.getActiveFilename;
    if ~isempty(activeFilename)
        candidateDirectories{end + 1} = fileparts(activeFilename);
    end
catch
    % MATLAB may be running without the desktop Editor.
end

scriptFullPath = mfilename('fullpath');
if ~isempty(scriptFullPath)
    candidateDirectories{end + 1} = fileparts(scriptFullPath);
end
resolvedScriptPath = which('plot_species_vs_eta');
if ~isempty(resolvedScriptPath)
    candidateDirectories{end + 1} = fileparts(resolvedScriptPath);
end
candidateDirectories{end + 1} = pwd;
candidateDirectories = unique(candidateDirectories, 'stable');

scriptDirectory = '';
for directoryIndex = 1:numel(candidateDirectories)
    candidateWorkbookPath = fullfile( ...
        candidateDirectories{directoryIndex}, inputFileName);
    candidateScriptPath = fullfile(candidateDirectories{directoryIndex}, ...
        'plot_species_vs_eta.m');
    if isfile(candidateWorkbookPath) && isfile(candidateScriptPath)
        scriptDirectory = candidateDirectories{directoryIndex};
        break
    end
end
if isempty(scriptDirectory)
    error('Could not locate %s beside plot_species_vs_eta.m.', ...
        inputFileName)
end
inputWorkbookPath = fullfile(scriptDirectory, inputFileName);

% Read the exact mean and standard-deviation columns used by the figure.
sheetName = 'Figure data';
dataTable = readtable(inputWorkbookPath, 'Sheet', sheetName, ...
    'VariableNamingRule', 'preserve');
variableNames = string(dataTable.Properties.VariableNames);
requiredVariables = ["eta", "mean_fraction_surviving", ...
    "standard_deviation"];
if any(~ismember(requiredVariables, variableNames))
    error(['Sheet %s must contain eta, mean_fraction_surviving, and ' ...
        'standard_deviation.'], sheetName)
end

etaValues = dataTable{:, variableNames == "eta"};
meanSurvivingFraction = dataTable{:, ...
    variableNames == "mean_fraction_surviving"};
stdSurvivingFraction = dataTable{:, ...
    variableNames == "standard_deviation"};

etaValues = etaValues(:);
meanSurvivingFraction = meanSurvivingFraction(:);
stdSurvivingFraction = stdSurvivingFraction(:);
if ~isnumeric(etaValues) || ~isnumeric(meanSurvivingFraction) || ...
        ~isnumeric(stdSurvivingFraction) || isempty(etaValues) || ...
        numel(meanSurvivingFraction) ~= numel(etaValues) || ...
        numel(stdSurvivingFraction) ~= numel(etaValues) || ...
        any(~isfinite(etaValues)) || ...
        any(~isfinite(meanSurvivingFraction)) || ...
        any(~isfinite(stdSurvivingFraction))
    error('The saved species-versus-eta data are invalid or incomplete.')
end
if any(diff(etaValues) < 0) || any(stdSurvivingFraction < 0)
    error('Eta must be nondecreasing and standard deviations nonnegative.')
end

% Match the original single-panel error-bar figure formatting.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 110;
markerSize = 30;
markerLineWidth = 3;
axesLineWidth = 1.5;
plotColor = lines(1);

survivingFractionMinimum = min( ...
    meanSurvivingFraction - stdSurvivingFraction, [], 'omitnan');
if ~isfinite(survivingFractionMinimum) || survivingFractionMinimum >= 1
    survivingFractionMinimum = 0.9;
end

figureHandle = figure(1);
clf(figureHandle)
set(figureHandle, 'Color', 'w', 'Position', [180 100 900 850], ...
    'Renderer', 'painters')
errorbar(etaValues, meanSurvivingFraction, stdSurvivingFraction, ...
    'Color', plotColor, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', markerSize, ...
    'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor, ...
    'LineWidth', markerLineWidth, 'CapSize', 10)
xlabel('$\eta$')
ylabel('Fraction of surviving species')
xlim([min(etaValues), max(etaValues)])
ylim([survivingFractionMinimum, 1])
xticks([0, 1, 2])
set(gca, 'FontSize', figureFontSize, 'LineWidth', axesLineWidth)
pbaspect([1 1 1])
box on

fprintf('Plotted %d eta values from %s.\n', ...
    numel(etaValues), inputWorkbookPath)
