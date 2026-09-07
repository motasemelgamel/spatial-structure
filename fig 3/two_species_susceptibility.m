% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Calculate the two-species steady states and susceptibility coefficient.

% Consumer-resource model parameters.
model.S = 2;
model.alpha = [5, 5];
model.k = [1, 1];
model.gammaR = 0.01;
model.gammaT = model.gammaR;
model.Rstar = 5;
model.N0 = 1;
model.R0 = model.Rstar;
model.T0 = model.Rstar;
model.extinctionThreshold = 1e-3;

% Parameter grids. Beta is logarithmically sampled over the established
% range, while theta retains the resolution of the two-species calculation.
betaValues = logspace(-5, 5, 100);
thetaValues = 0.05:0.01:0.95;

% Simulation settings follow the adaptive explicit-Euler safety convention.
tMax = 3000;
explicitSafetyFactor = 0.5;
maximumDtValues = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* betaValues));
minimumSteadyTime = 1;
steadyCheckInterval = 1;
steadyTolerance = 1e-8;

nTheta = numel(thetaValues);
nBeta = numel(betaValues);
N1 = zeros(nTheta, nBeta);
N2 = zeros(nTheta, nBeta);
T1 = zeros(nTheta, nBeta);
T2 = zeros(nTheta, nBeta);
R = zeros(nTheta, nBeta);
steadyStateConfirmed = false(nBeta, 1);
minimumDtUsed = zeros(nBeta, 1);
maximumDtUsed = zeros(nBeta, 1);
numSteps = zeros(nBeta, 1);

% Every beta value is independent, so calculate them concurrently.
parfor betaIndex = 1:nBeta
    [currentN, currentR, currentT, currentSteadyStateConfirmed, ...
        currentMinimumDt, currentMaximumDt, currentNumSteps] = ...
        simulateToSteadyState(betaValues(betaIndex), thetaValues, tMax, ...
        maximumDtValues(betaIndex), explicitSafetyFactor, ...
        minimumSteadyTime, steadyCheckInterval, steadyTolerance, model);

    N1(:, betaIndex) = currentN(:, 1);
    N2(:, betaIndex) = currentN(:, 2);
    T1(:, betaIndex) = currentT(:, 1);
    T2(:, betaIndex) = currentT(:, 2);
    R(:, betaIndex) = currentR;
    steadyStateConfirmed(betaIndex) = currentSteadyStateConfirmed;
    minimumDtUsed(betaIndex) = currentMinimumDt;
    maximumDtUsed(betaIndex) = currentMaximumDt;
    numSteps(betaIndex) = currentNumSteps;
end

% Report any parameter values that did not satisfy the steady-state test.
unconfirmedIndices = find(~steadyStateConfirmed);
if ~isempty(unconfirmedIndices)
    warning(['A steady state was not confirmed before tMax for %d beta ' ...
        'values: %s'], numel(unconfirmedIndices), ...
        mat2str(betaValues(unconfirmedIndices), 4))
end

fprintf('Completed %d beta values over [%.1e, %.1e].\n', ...
    nBeta, min(betaValues), max(betaValues))
fprintf('Overall timestep range: [%.6g, %.6g].\n', ...
    min(minimumDtUsed), max(maximumDtUsed))
fprintf('Total adaptive Euler steps: %d.\n', sum(numSteps))

% Calculate only the susceptibility coefficient used in the output figure.
susceptibility11 = calculateSusceptibility11( ...
    N1, N2, T1, T2, R, thetaValues, model.alpha, betaValues, ...
    model.gammaR, model.gammaT);

% Prepare finite heatmap values without changing the calculated results.
susceptibilityMagnitude = abs(susceptibility11);
validSusceptibility = isfinite(susceptibilityMagnitude) & ...
    susceptibilityMagnitude > 0;
if ~any(validSusceptibility, 'all')
    error('No finite positive susceptibility values are available to plot.')
end
susceptibilityReference = min( ...
    susceptibilityMagnitude(validSusceptibility));
susceptibilityForPlot = nan(size(susceptibilityMagnitude));
susceptibilityForPlot(validSusceptibility) = log10( ...
    susceptibilityMagnitude(validSusceptibility) ./ ...
    susceptibilityReference);

% Display extinct populations at the extinction floor on the log scale.
logN1ForPlot = log10(max(N1, model.extinctionThreshold));
logN2ForPlot = log10(max(N2, model.extinctionThreshold));
logBetaValues = log10(betaValues);
abundanceColorLimits = [min([logN1ForPlot(:); logN2ForPlot(:)]), ...
    max([logN1ForPlot(:); logN2ForPlot(:)])];
if abundanceColorLimits(1) == abundanceColorLimits(2)
    abundanceColorLimits = abundanceColorLimits + [-0.5, 0.5];
end

% Apply the shared figure style used in the other parameter-sweep figures.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
figureFontSize = 60;
titleFontSize = 46;
colorbarFontSize = 42;
axesLineWidth = 1.5;

% Figure 1: normalized susceptibility coefficient.
figure(1)
clf
set(gcf, 'Color', 'w', 'Position', [180 100 900 850])
susceptibilityAxes = axes;
formatHeatmapAxes(susceptibilityAxes, thetaValues, logBetaValues, ...
    susceptibilityForPlot.', figureFontSize, axesLineWidth)
title(susceptibilityAxes, ...
    '$\log_{10}(|\chi_{11}|/\min |\chi_{11}|)$', ...
    'FontSize', titleFontSize)
formatHeatmapColorbar(susceptibilityAxes, colorbarFontSize, axesLineWidth)

% Figure 2: equilibrium abundance of each species.
figure(2)
clf
set(gcf, 'Color', 'w', 'Position', [50 80 1900 900])
abundanceLayout = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

firstAbundanceAxes = nexttile(abundanceLayout, 1);
formatHeatmapAxes(firstAbundanceAxes, thetaValues, logBetaValues, ...
    logN1ForPlot.', figureFontSize, axesLineWidth)
clim(firstAbundanceAxes, abundanceColorLimits)
title(firstAbundanceAxes, '$\log_{10}(N_1)$', ...
    'FontSize', titleFontSize)

secondAbundanceAxes = nexttile(abundanceLayout, 2);
formatHeatmapAxes(secondAbundanceAxes, thetaValues, logBetaValues, ...
    logN2ForPlot.', figureFontSize, axesLineWidth)
clim(secondAbundanceAxes, abundanceColorLimits)
title(secondAbundanceAxes, '$\log_{10}(N_2)$', ...
    'FontSize', titleFontSize)
abundanceColorbar = formatHeatmapColorbar( ...
    secondAbundanceAxes, colorbarFontSize, axesLineWidth);
abundanceColorbar.Layout.Tile = 'east';

function [N, R, T, steadyStateConfirmed, minimumDtUsed, ...
    maximumDtUsed, numSteps] = simulateToSteadyState( ...
    beta, thetaValues, tMax, maximumDt, explicitSafetyFactor, ...
    minimumSteadyTime, steadyCheckInterval, steadyTolerance, model)

% Simulate all theta values for one beta using a shared adaptive timestep.
thetaColumn = thetaValues(:);
nTheta = numel(thetaColumn);
N = model.N0 .* ones(nTheta, model.S);
R = model.R0 .* ones(nTheta, 1);
T = model.T0 .* ones(nTheta, model.S);

currentTime = 0;
nextSteadyCheck = minimumSteadyTime;
steadyStateConfirmed = false;
minimumDtUsed = inf;
maximumDtUsed = 0;
numSteps = 0;

while currentTime < tMax
    [dN, dR, dT, perCapitaGrowth, uptake] = calculateRates( ...
        N, R, T, beta, thetaColumn, model);

    % Test the maximum normalized derivative at fixed physical intervals.
    if currentTime >= nextSteadyCheck
        maximumRelativeRate = calculateMaximumRelativeRate( ...
            N, R, T, dN, dR, dT);
        if maximumRelativeRate < steadyTolerance
            steadyStateConfirmed = true;
            break
        end
        nextSteadyCheck = currentTime + steadyCheckInterval;
    end

    % Bound the explicit timestep using exchange, consumption, and growth.
    localRemovalRate = beta + model.gammaT .* ...
        thetaColumn .* uptake;
    sharedRemovalRate = 1 + model.S .* beta + model.gammaR .* ...
        (1 - thetaColumn) .* sum(uptake, 2);
    survivingNow = N > model.extinctionThreshold;
    if any(survivingNow, 'all')
        demographicRate = max( ...
            abs(perCapitaGrowth(survivingNow)), [], 'all');
    else
        demographicRate = 0;
    end

    maximumCurrentRate = max([localRemovalRate(:); ...
        sharedRemovalRate(:); demographicRate]);
    if ~isfinite(maximumCurrentRate) || maximumCurrentRate <= 0
        error(['An invalid model rate occurred at t = %.6g for ' ...
            'beta = %.6g.'], currentTime, beta)
    end

    dt = min(maximumDt, ...
        explicitSafetyFactor ./ maximumCurrentRate);
    dt = min(dt, tMax - currentTime);
    if ~isfinite(dt) || dt <= 0
        error(['An invalid timestep occurred at t = %.6g for ' ...
            'beta = %.6g.'], currentTime, beta)
    end

    % Advance the consumer-resource model explicitly.
    NNew = N + dt .* dN;
    RNew = R + dt .* dR;
    TNew = T + dt .* dT;

    % Force extinct populations to zero immediately and permanently.
    NNew(NNew < model.extinctionThreshold) = 0;

    % Reject invalid numerical states instead of plotting them silently.
    if any(~isfinite(NNew), 'all') || any(~isfinite(RNew), 'all') || ...
            any(~isfinite(TNew), 'all')
        error(['A nonfinite state occurred at t = %.6g for ' ...
            'beta = %.6g.'], currentTime + dt, beta)
    end
    if any(RNew < 0, 'all') || any(TNew < 0, 'all')
        error(['A negative resource state occurred at t = %.6g for ' ...
            'beta = %.6g.'], currentTime + dt, beta)
    end

    % Commit the completed adaptive Euler step.
    N = NNew;
    R = RNew;
    T = TNew;
    currentTime = currentTime + dt;
    numSteps = numSteps + 1;
    minimumDtUsed = min(minimumDtUsed, dt);
    maximumDtUsed = max(maximumDtUsed, dt);
end

% Evaluate the final state when the trajectory reaches tMax between checks.
if ~steadyStateConfirmed
    [dN, dR, dT] = calculateRates( ...
        N, R, T, beta, thetaColumn, model);
    maximumRelativeRate = calculateMaximumRelativeRate( ...
        N, R, T, dN, dR, dT);
    steadyStateConfirmed = maximumRelativeRate < steadyTolerance;
end
end

function [dN, dR, dT, perCapitaGrowth, uptake] = ...
    calculateRates(N, R, T, beta, thetaColumn, model)

% Calculate the deterministic consumer-resource rates for all theta values.
oneMinusTheta = 1 - thetaColumn;
resourceUse = thetaColumn .* T + oneMinusTheta .* R;
perCapitaGrowth = model.alpha .* resourceUse - model.k;
uptake = model.alpha .* N;

dN = perCapitaGrowth .* N;
dR = model.Rstar - R - beta .* ...
    (model.S .* R - sum(T, 2)) - model.gammaR .* R .* ...
    oneMinusTheta .* sum(uptake, 2);
dT = beta .* (R - T) - model.gammaT .* ...
    thetaColumn .* T .* uptake;
end

function maximumRelativeRate = calculateMaximumRelativeRate( ...
    N, R, T, dN, dR, dT)

% Normalize derivatives by the magnitude of their corresponding state.
relativeRates = [ ...
    abs(dN(:)) ./ max(1, abs(N(:))); ...
    abs(dR(:)) ./ max(1, abs(R(:))); ...
    abs(dT(:)) ./ max(1, abs(T(:)))];
maximumRelativeRate = max(relativeRates);
end

function susceptibility11 = calculateSusceptibility11( ...
    N1, N2, T1, T2, R, thetaValues, alpha, betaValues, ...
    gammaR, gammaT)

% Calculate the first diagonal coefficient of the susceptibility matrix.
thetaColumn = thetaValues(:);
betaRow = betaValues(:).';
oneMinusTheta = 1 - thetaColumn;

firstDenominator = betaRow + ...
    alpha(1) .* gammaT .* thetaColumn .* N1;
secondDenominator = betaRow + ...
    alpha(2) .* gammaT .* thetaColumn .* N2;

A1 = alpha(1) .* oneMinusTheta + ...
    alpha(1) .* thetaColumn .* betaRow ./ firstDenominator;
A2 = -(alpha(1).^2 .* thetaColumn.^2 .* gammaT .* T1) ./ ...
    firstDenominator;
B1 = alpha(2) .* oneMinusTheta + ...
    alpha(2) .* thetaColumn .* betaRow ./ secondDenominator;
B2 = -(alpha(2).^2 .* thetaColumn.^2 .* gammaT .* T2) ./ ...
    secondDenominator;

C1 = -1 - 2 .* betaRow + betaRow.^2 ./ firstDenominator + ...
    betaRow.^2 ./ secondDenominator - gammaR .* oneMinusTheta .* ...
    (alpha(1) .* N1 + alpha(2) .* N2);
C2 = -(alpha(1) .* gammaT .* thetaColumn .* T1 .* betaRow) ./ ...
    firstDenominator - gammaR .* oneMinusTheta .* alpha(1) .* R;
C3 = -(alpha(2) .* gammaT .* thetaColumn .* T2 .* betaRow) ./ ...
    secondDenominator - gammaR .* oneMinusTheta .* alpha(2) .* R;

susceptibilityDenominator = A2 .* B2 .* C1 - ...
    A1 .* B2 .* C2 - A2 .* B1 .* C3;
susceptibility11 = (B2 .* C1 - B1 .* C3) ./ ...
    susceptibilityDenominator;
end

function formatHeatmapAxes(axesHandle, thetaValues, logBetaValues, ...
    heatmapValues, figureFontSize, axesLineWidth)

% Plot uniformly sampled log10(beta) coordinates without axis distortion.
imagesc(axesHandle, thetaValues, logBetaValues, heatmapValues)
set(axesHandle, 'YDir', 'normal', 'FontSize', figureFontSize, ...
    'LineWidth', axesLineWidth, 'TickDir', 'out', 'Layer', 'top')
xlabel(axesHandle, '$\theta$')
ylabel(axesHandle, '$\beta$')
xlim(axesHandle, [min(thetaValues), max(thetaValues)])
ylim(axesHandle, [min(logBetaValues), max(logBetaValues)])
xticks(axesHandle, [0.05, 0.5, 0.95])

minimumLogBeta = floor(min(logBetaValues));
maximumLogBeta = ceil(max(logBetaValues));
logBetaTicks = minimumLogBeta:2:maximumLogBeta;
if logBetaTicks(end) < maximumLogBeta
    logBetaTicks(end + 1) = maximumLogBeta;
end
yticks(axesHandle, logBetaTicks)
yticklabels(axesHandle, compose('$10^{%d}$', logBetaTicks))
colormap(axesHandle, parula(256))
pbaspect(axesHandle, [1 1 1])
box(axesHandle, 'on')
end

function colorbarHandle = formatHeatmapColorbar( ...
    axesHandle, colorbarFontSize, axesLineWidth)

% Apply the established typography and line weight to one colorbar.
colorbarHandle = colorbar(axesHandle);
colorbarHandle.FontSize = colorbarFontSize;
colorbarHandle.LineWidth = axesLineWidth;
colorbarHandle.TickLabelInterpreter = 'latex';
end
