% Disclaimer: This script was written and revised with assistance from
% OpenAI Codex. The author remains responsible for verifying all results.

clearvars
close all
clc

% Calculate the community-averaged deterministic surviving fraction across
% adaptation rates.

% Set the reproducible master seed used to sample independent communities.
masterSeed = 1;
traitStream = RandStream('Threefry', 'Seed', masterSeed);

% Consumer-resource model parameters.
model.S = 50;
model.k = 1;
model.gamma = 0.01;
model.Rstar = 5;
model.N0 = 1;
model.R0 = model.Rstar;
model.T0 = model.Rstar;
model.extinctionThreshold = 1e-3;

% Fixed diffusion rate and the uniformly sampled adaptation-rate grid.
beta = 1;
etaValues = linspace(0, 2, 100);

% Sample 20 independent communities once and reuse the same communities at
% every eta value. This paired design isolates the effect of eta from
% differences in the sampled alpha values.
nReplicates = 20;
alphaByReplicate = zeros(nReplicates, model.S);
for replicateIndex = 1:nReplicates
    traitStream.Substream = replicateIndex;
    alphaByReplicate(replicateIndex, :) = ...
        1 + 4 .* rand(traitStream, 1, model.S);
end

% Simulation settings retained from the adaptive dynamics calculation.
tMax = 3000;
explicitSafetyFactor = 0.5;
maximumDt = min(0.001, explicitSafetyFactor ./ ...
    ((model.S + 1) .* beta));

nEta = numel(etaValues);
numberOfSimulations = nEta .* nReplicates;
survivingFractionByReplicate = zeros(numberOfSimulations, 1);
minimumDtUsed = zeros(numberOfSimulations, 1);
maximumDtUsed = zeros(numberOfSimulations, 1);
numSteps = zeros(numberOfSimulations, 1);

% Simulate every eta-community pair concurrently across available workers.
parfor simulationIndex = 1:numberOfSimulations
    [etaIndex, replicateIndex] = ind2sub( ...
        [nEta, nReplicates], simulationIndex);
    eta = etaValues(etaIndex);
    replicateModel = model;
    replicateModel.alpha = alphaByReplicate(replicateIndex, :);

    [survivingFractionByReplicate(simulationIndex), ...
        minimumDtUsed(simulationIndex), ...
        maximumDtUsed(simulationIndex), numSteps(simulationIndex)] = ...
        simulateFinalSurvivingFraction( ...
        eta, beta, tMax, maximumDt, explicitSafetyFactor, replicateModel);

    fprintf(['Completed simulation %d of %d: eta = %.6g, ' ...
        'community = %d, fraction = %.6g.\n'], ...
        simulationIndex, numberOfSimulations, eta, replicateIndex, ...
        survivingFractionByReplicate(simulationIndex))
end

% Recover the eta-by-community matrix and calculate ensemble statistics.
survivingFractionByReplicate = reshape( ...
    survivingFractionByReplicate, [nEta, nReplicates]);
meanSurvivingFraction = mean( ...
    survivingFractionByReplicate, 2, 'omitnan');
stdSurvivingFraction = std( ...
    survivingFractionByReplicate, 0, 2, 'omitnan');

% Match the single-panel figure formatting used for the parameter sweeps.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figureFontSize = 60;
markerSize = 15;
markerLineWidth = 3;
axesLineWidth = 1.5;
plotColor = lines(1);

% Include the lower standard-deviation bars while retaining a visible range
% if every replicate has complete survival throughout the sweep.
survivingFractionMinimum = min( ...
    meanSurvivingFraction - stdSurvivingFraction, [], 'omitnan');
if ~isfinite(survivingFractionMinimum) || survivingFractionMinimum >= 1
    survivingFractionMinimum = 0.9;
end

% The only output figure: community mean with one-standard-deviation bars.
figure(1)
clf
set(gcf, 'Color', 'w', 'Position', [180 100 900 850])
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

% Report completion and timestep diagnostics across all simulations.
fprintf(['All %d simulations completed for %d eta values and %d ' ...
    'communities.\n'], numberOfSimulations, nEta, nReplicates)
fprintf('Overall timestep range: [%.6g, %.6g].\n', ...
    min(minimumDtUsed), max(maximumDtUsed))
fprintf('Total adaptive Euler steps: %d.\n', sum(numSteps))

function [survivingFraction, minimumDtUsed, maximumDtUsed, numSteps] = ...
    simulateFinalSurvivingFraction( ...
    eta, beta, tMax, maximumDt, explicitSafetyFactor, model)

% Simulate one deterministic eta value and retain only its final richness.
N = model.N0 .* ones(1, model.S);
R = model.R0;
T = model.T0 .* ones(1, model.S);
theta = zeros(1, model.S);

currentTime = 0;
minimumDtUsed = inf;
maximumDtUsed = 0;
numSteps = 0;

while currentTime < tMax
    % Calculate each species' current per-capita growth rate.
    resourceUse = theta .* T + (1 - theta) .* R;
    perCapitaGrowth = model.alpha .* resourceUse - model.k;

    % Calculate the deterministic consumer-resource rates.
    dN = perCapitaGrowth .* N;
    dT = beta .* (R - T) - ...
        model.alpha .* model.gamma .* theta .* T .* N;
    dR = model.Rstar - R - beta .* ...
        (model.S .* R - sum(T)) - model.gamma .* R .* ...
        sum(model.alpha .* (1 - theta) .* N);

    % Bound the explicit timestep using the current resource, demographic,
    % and adaptation rates. Extinct species have frozen N_i and theta_i.
    survivingNow = N > model.extinctionThreshold;
    localRemovalRate = beta + ...
        model.alpha .* model.gamma .* theta .* N;
    sharedRemovalRate = 1 + model.S .* beta + ...
        model.gamma .* sum(model.alpha .* (1 - theta) .* N);

    if any(survivingNow)
        demographicRate = max(abs(perCapitaGrowth(survivingNow)));
        adaptationRate = max(eta .* ...
            abs(perCapitaGrowth(survivingNow)));
    else
        demographicRate = 0;
        adaptationRate = 0;
    end

    maximumCurrentRate = max([localRemovalRate, sharedRemovalRate, ...
        demographicRate, adaptationRate]);

    if ~isfinite(maximumCurrentRate) || maximumCurrentRate <= 0
        error(['An invalid adaptive rate occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end

    dt = min(maximumDt, ...
        explicitSafetyFactor ./ maximumCurrentRate);
    dt = min(dt, tMax - currentTime);

    if ~isfinite(dt) || dt <= 0
        error(['An invalid timestep occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end

    % Advance abundances and resources explicitly.
    NNew = N + dt .* dN;
    TNew = T + dt .* dT;
    RNew = R + dt .* dR;

    % Force populations below the extinction threshold permanently to zero.
    NNew(NNew < model.extinctionThreshold) = 0;

    % Advance theta semi-implicitly with growth held fixed during the step.
    switchLocal = eta .* max(-perCapitaGrowth, 0);
    switchGlobal = eta .* max(perCapitaGrowth, 0);
    thetaNew = theta;

    negativeGrowth = perCapitaGrowth < 0 & NNew > 0;
    positiveGrowth = perCapitaGrowth > 0 & NNew > 0;

    thetaNew(negativeGrowth) = ...
        (theta(negativeGrowth) + dt .* switchLocal(negativeGrowth)) ...
        ./ (1 + dt .* switchLocal(negativeGrowth));
    thetaNew(positiveGrowth) = theta(positiveGrowth) ...
        ./ (1 + dt .* switchGlobal(positiveGrowth));

    % Extinct species retain their final theta, and roundoff is bounded.
    thetaNew(NNew == 0) = theta(NNew == 0);
    thetaNew = min(max(thetaNew, 0), 1);

    % Commit the completed adaptive step.
    N = NNew;
    T = TNew;
    R = RNew;
    theta = thetaNew;
    currentTime = currentTime + dt;
    numSteps = numSteps + 1;
    minimumDtUsed = min(minimumDtUsed, dt);
    maximumDtUsed = max(maximumDtUsed, dt);

    % Stop instead of retaining a numerically invalid trajectory.
    if any(~isfinite(N), 'all') || ~isfinite(R) || ...
            any(~isfinite(T), 'all') || any(~isfinite(theta), 'all')
        error(['A nonfinite state occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end
    if R < 0 || any(T < 0, 'all')
        error(['A negative resource state occurred at t = %.6g ' ...
            'for eta = %.6g.'], currentTime, eta)
    end
end

survivingFraction = nnz(N > model.extinctionThreshold) ./ model.S;
end
