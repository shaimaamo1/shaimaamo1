%% CML Numerical Analysis Project
% MATLAB-Based Numerical Simulation of BCR-ABL1 Expression Dynamics
% and Targeted Therapy Response in Chronic Myeloid Leukemia
%
% This script includes:
% 1. Gene expression data loading
% 2. Basic preprocessing
% 3. Differential expression analysis
% 4. Scientific visualization
% 5. Numerical simulation of CML dynamics
% 6. Euler, RK4, and ode45 comparison
% 7. Treatment-response scenarios
%
% Place this file in the same folder as the dataset, then press Run.
% If no dataset is found, a small demo dataset is generated automatically
% so that the full workflow can still run without errors.

clear;
clc;
close all;

%% =========================
% USER SETTINGS
%% =========================

outputFolder = "Project_Output";
figFolder = fullfile(outputFolder, "Figures");

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

possibleFiles = [
    "GSE33075_expression.csv"
    "GSE33075_expression.xlsx"
    "gene_expression.csv"
    "gene_expression.xlsx"
    "expression_data.csv"
    "expression_data.xlsx"
];

groupA = "CML";
groupB = "Normal";

logFC_cutoff = 1;
p_cutoff = 0.05;
numTopGenes = 10;

%% =========================
% STEP 1: LOAD EXPRESSION DATA
%% =========================

fprintf('\nSTEP 1: Loading expression data...\n');
[dataTable, usedFile, isDemo] = loadExpressionData(possibleFiles);

if isDemo
    fprintf('No dataset was found. Demo data was generated for testing.\n');
else
    fprintf('Dataset loaded from: %s\n', usedFile);
end

varNames = string(dataTable.Properties.VariableNames);
possibleGeneCols = ["Gene", "GeneSymbol", "Gene_Symbol", "Symbol", "gene", "ID", "ProbeID"];
geneColIndex = find(ismember(varNames, possibleGeneCols), 1);

if isempty(geneColIndex)
    geneColIndex = 1;
end

geneNames = string(dataTable{:, geneColIndex});
exprTable = dataTable;
exprTable(:, geneColIndex) = [];

sampleNames = string(exprTable.Properties.VariableNames);
expressionMatrix = double(table2array(exprTable));

fprintf('Genes loaded: %d\n', size(expressionMatrix, 1));
fprintf('Samples loaded: %d\n', size(expressionMatrix, 2));

%% =========================
% STEP 2: CREATE METADATA
%% =========================

fprintf('\nSTEP 2: Creating metadata...\n');
metadata = createMetadata(sampleNames, groupA, groupB);
disp(metadata);

%% =========================
% STEP 3: DATA CLEANING
%% =========================

fprintf('\nSTEP 3: Cleaning data...\n');

validGeneNames = geneNames ~= "" & ~ismissing(geneNames);
geneNames = geneNames(validGeneNames);
expressionMatrix = expressionMatrix(validGeneNames, :);

expressionMatrix = fillMissingByGeneMean(expressionMatrix);

maxValue = max(expressionMatrix(:), [], 'omitnan');
if maxValue > 50
    expressionMatrix = log2(expressionMatrix + 1);
    fprintf('Log2 transformation was applied.\n');
else
    fprintf('Data already appears log-transformed.\n');
end

variancePerGene = var(expressionMatrix, 0, 2, 'omitnan');
keepGenes = variancePerGene > 0;
expressionMatrix = expressionMatrix(keepGenes, :);
geneNames = geneNames(keepGenes);

fprintf('Genes after cleaning: %d\n', length(geneNames));

%% =========================
% STEP 4: QUALITY CONTROL FIGURES
%% =========================

fprintf('\nSTEP 4: Creating QC figures...\n');

figure('Name', 'Sample Expression Boxplot');
boxplot(expressionMatrix, 'Labels', sampleNames);
ylabel('Expression value');
title('Expression Distribution Across Samples');
xtickangle(45);
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_1_sample_expression_boxplot.png'));

sampleMeans = mean(expressionMatrix, 1, 'omitnan');
sampleStd = std(expressionMatrix, 0, 1, 'omitnan');

figure('Name', 'Mean Expression per Sample');
bar(sampleMeans);
hold on;
errorbar(1:length(sampleMeans), sampleMeans, sampleStd, 'k.', 'LineWidth', 1);
hold off;
set(gca, 'XTick', 1:length(sampleNames), 'XTickLabel', sampleNames);
xtickangle(45);
ylabel('Mean expression');
title('Mean Expression per Sample');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_2_sample_mean_expression.png'));

%% =========================
% STEP 5: DIFFERENTIAL EXPRESSION ANALYSIS
%% =========================

fprintf('\nSTEP 5: Differential expression analysis...\n');
fprintf('Comparison: %s vs %s\n', groupA, groupB);

idxA = metadata.Group == groupA;
idxB = metadata.Group == groupB;

if sum(idxA) < 2 || sum(idxB) < 2
    warning('Group labels were not clear. Samples were split automatically.');
    halfPoint = floor(size(expressionMatrix, 2) / 2);
    idxB = false(1, size(expressionMatrix, 2));
    idxA = false(1, size(expressionMatrix, 2));
    idxB(1:halfPoint) = true;
    idxA(halfPoint+1:end) = true;
end

meanA = mean(expressionMatrix(:, idxA), 2, 'omitnan');
meanB = mean(expressionMatrix(:, idxB), 2, 'omitnan');
logFC = meanA - meanB;
pValues = zeros(size(expressionMatrix, 1), 1);

for i = 1:size(expressionMatrix, 1)
    valuesA = expressionMatrix(i, idxA);
    valuesB = expressionMatrix(i, idxB);
    [~, p] = ttest2(valuesA, valuesB, 'Vartype', 'unequal');
    pValues(i) = p;
end

adjPValues = benjaminiHochberg(pValues);

regulation = strings(length(geneNames), 1);
regulation(:) = "Not significant";
regulation(logFC >= logFC_cutoff & adjPValues < p_cutoff) = "Upregulated";
regulation(logFC <= -logFC_cutoff & adjPValues < p_cutoff) = "Downregulated";

DEG_Table = table(geneNames, meanA, meanB, logFC, pValues, adjPValues, regulation, ...
    'VariableNames', {'Gene', 'Mean_GroupA', 'Mean_GroupB', 'logFC', 'PValue', 'AdjustedPValue', 'Regulation'});

DEG_Table = sortrows(DEG_Table, 'AdjustedPValue', 'ascend');
writetable(DEG_Table, fullfile(outputFolder, 'Differential_Expression_Results.xlsx'));

numUp = sum(DEG_Table.Regulation == "Upregulated");
numDown = sum(DEG_Table.Regulation == "Downregulated");

fprintf('Upregulated genes: %d\n', numUp);
fprintf('Downregulated genes: %d\n', numDown);

%% =========================
% STEP 6: DEG FIGURES
%% =========================

fprintf('\nSTEP 6: Creating DEG figures...\n');

negLogP = -log10(adjPValues);
notSig = regulation == "Not significant";
upReg = regulation == "Upregulated";
downReg = regulation == "Downregulated";

figure('Name', 'Volcano Plot');
scatter(logFC(notSig), negLogP(notSig), 18, 'filled');
hold on;
scatter(logFC(upReg), negLogP(upReg), 22, 'filled');
scatter(logFC(downReg), negLogP(downReg), 22, 'filled');
xline(logFC_cutoff, '--');
xline(-logFC_cutoff, '--');
yline(-log10(p_cutoff), '--');
hold off;
xlabel('log2 Fold Change');
ylabel('-log10 Adjusted P-value');
title('Volcano Plot of Differentially Expressed Genes');
legend({'Not significant', 'Upregulated', 'Downregulated'}, 'Location', 'best');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_3_volcano_plot.png'));

sigGenes = DEG_Table(DEG_Table.Regulation ~= "Not significant", :);
if height(sigGenes) > 0
    topN = min(numTopGenes, height(sigGenes));
    topGenes = sigGenes(1:topN, :);
else
    topN = min(numTopGenes, height(DEG_Table));
    topGenes = DEG_Table(1:topN, :);
end

figure('Name', 'Top Differentially Expressed Genes');
bar(topGenes.logFC);
set(gca, 'XTick', 1:height(topGenes), 'XTickLabel', topGenes.Gene);
xtickangle(45);
ylabel('log2 Fold Change');
title('Top Differentially Expressed Genes');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_4_top_DE_genes_barplot.png'));

[~, geneIndexForHeatmap] = ismember(topGenes.Gene, geneNames);
heatmapData = expressionMatrix(geneIndexForHeatmap, :);
heatmapData = zscore(heatmapData, 0, 2);

figure('Name', 'Top Genes Heatmap');
imagesc(heatmapData);
colorbar;
set(gca, 'XTick', 1:length(sampleNames), 'XTickLabel', sampleNames);
set(gca, 'YTick', 1:height(topGenes), 'YTickLabel', topGenes.Gene);
xtickangle(45);
title('Heatmap of Top Differentially Expressed Genes');
xlabel('Samples');
ylabel('Genes');
saveas(gcf, fullfile(figFolder, 'Figure_5_top_genes_heatmap.png'));

%% =========================
% STEP 7: BIOLOGICAL SUMMARY FILE
%% =========================

fprintf('\nSTEP 7: Writing project summary...\n');

summaryFile = fullfile(outputFolder, 'Project_Summary.txt');
fid = fopen(summaryFile, 'w');

fprintf(fid, 'Numerical Analysis Project Summary\n');
fprintf(fid, '==================================\n\n');
fprintf(fid, 'Comparison: %s vs %s\n', groupA, groupB);
fprintf(fid, 'Genes analyzed: %d\n', length(geneNames));
fprintf(fid, 'Samples analyzed: %d\n', length(sampleNames));
fprintf(fid, 'Upregulated genes: %d\n', numUp);
fprintf(fid, 'Downregulated genes: %d\n\n', numDown);
fprintf(fid, 'Top genes by adjusted p-value:\n');

for i = 1:min(numTopGenes, height(DEG_Table))
    fprintf(fid, '%d. %s | logFC = %.3f | adjusted p-value = %.4g | %s\n', ...
        i, DEG_Table.Gene(i), DEG_Table.logFC(i), DEG_Table.AdjustedPValue(i), DEG_Table.Regulation(i));
end

fclose(fid);

%% =========================
% STEP 8: CML MATHEMATICAL MODEL
%% =========================

fprintf('\nSTEP 8: Running CML mathematical model...\n');

params.rL = 0.42;      % leukemic cell growth rate
params.KL = 100;       % leukemic carrying capacity
params.dL = 0.04;      % leukemic natural death rate
params.rH = 0.30;      % healthy cell growth rate
params.KH = 100;       % healthy carrying capacity
params.alpha = 0.004;  % competition effect of leukemic cells
params.beta = 0.18;    % mRNA production from leukemic cells
params.gamma = 0.25;   % mRNA degradation rate
params.u = 0.35;       % treatment effect strength

Y0 = [15; 85; 8];      % [leukemic cells; healthy cells; mRNA signal]
t0 = 0;
tEnd = 60;
h = 0.5;
timeVector = t0:h:tEnd;

[tEuler, YEuler] = eulerMethod(@(t, Y) cmlModel(t, Y, params), timeVector, Y0);
[tRK4, YRK4] = rk4Method(@(t, Y) cmlModel(t, Y, params), timeVector, Y0);
[tOde45, YOde45] = ode45(@(t, Y) cmlModel(t, Y, params), [t0 tEnd], Y0);

%% =========================
% STEP 9: MODEL FIGURES
%% =========================

fprintf('\nSTEP 9: Creating model figures...\n');

figure('Name', 'CML Dynamics using RK4');
plot(tRK4, YRK4(:, 1), 'LineWidth', 2);
hold on;
plot(tRK4, YRK4(:, 2), 'LineWidth', 2);
plot(tRK4, YRK4(:, 3), 'LineWidth', 2);
hold off;
xlabel('Time');
ylabel('Model value');
title('CML Dynamics Under Treatment using RK4');
legend({'Leukemic cells', 'Healthy cells', 'mRNA expression'}, 'Location', 'best');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_6_CML_dynamics_RK4.png'));

figure('Name', 'Numerical Method Comparison');
plot(tEuler, YEuler(:, 1), '--', 'LineWidth', 1.5);
hold on;
plot(tRK4, YRK4(:, 1), 'LineWidth', 2);
plot(tOde45, YOde45(:, 1), ':', 'LineWidth', 2);
hold off;
xlabel('Time');
ylabel('Leukemic cell population');
title('Numerical Method Comparison for Leukemic Cells');
legend({'Euler', 'RK4', 'ode45'}, 'Location', 'best');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_7_numerical_methods_comparison.png'));

%% =========================
% STEP 10: TREATMENT SCENARIOS
%% =========================

fprintf('\nSTEP 10: Comparing treatment scenarios...\n');

treatmentValues = [0, 0.15, 0.35, 0.55];
figure('Name', 'Treatment Scenario Comparison');
hold on;
scenarioSummary = table();

for k = 1:length(treatmentValues)
    tempParams = params;
    tempParams.u = treatmentValues(k);
    [tScenario, YScenario] = rk4Method(@(t, Y) cmlModel(t, Y, tempParams), timeVector, Y0);
    plot(tScenario, YScenario(:, 1), 'LineWidth', 2);

    scenarioSummary.Treatment(k, 1) = treatmentValues(k);
    scenarioSummary.Final_Leukemic_Cells(k, 1) = YScenario(end, 1);
    scenarioSummary.Final_Healthy_Cells(k, 1) = YScenario(end, 2);
    scenarioSummary.Final_mRNA(k, 1) = YScenario(end, 3);
end

hold off;
xlabel('Time');
ylabel('Leukemic cell population');
title('Effect of Treatment Strength on Leukemic Cell Dynamics');
legend({'u = 0', 'u = 0.15', 'u = 0.35', 'u = 0.55'}, 'Location', 'best');
grid on;
saveas(gcf, fullfile(figFolder, 'Figure_8_treatment_scenarios.png'));

writetable(scenarioSummary, fullfile(outputFolder, 'Treatment_Scenario_Summary.xlsx'));

fprintf('\nPROJECT FINISHED SUCCESSFULLY.\n');
fprintf('All outputs were saved in: %s\n', outputFolder);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [dataTable, usedFile, isDemo] = loadExpressionData(possibleFiles)
    usedFile = "";
    isDemo = false;

    for i = 1:length(possibleFiles)
        currentFile = possibleFiles(i);
        if exist(currentFile, 'file')
            usedFile = currentFile;
            [~, ~, ext] = fileparts(currentFile);

            if strcmpi(ext, '.csv') || strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
                dataTable = readtable(currentFile, 'VariableNamingRule', 'preserve');
                dataTable.Properties.VariableNames = matlab.lang.makeValidName(dataTable.Properties.VariableNames);
                return;
            end
        end
    end

    isDemo = true;
    rng(7);
    nGenes = 120;
    nNormal = 4;
    nCML = 4;

    geneNames = strings(nGenes, 1);
    for g = 1:nGenes
        geneNames(g) = "Gene_" + g;
    end

    normalData = 8 + randn(nGenes, nNormal);
    cmlData = 8 + randn(nGenes, nCML);
    cmlData(1:12, :) = cmlData(1:12, :) + 2.2;
    cmlData(13:24, :) = cmlData(13:24, :) - 2.0;

    expressionData = [normalData, cmlData];
    sampleNames = ["Normal_1", "Normal_2", "Normal_3", "Normal_4", ...
                   "CML_1", "CML_2", "CML_3", "CML_4"];

    dataTable = array2table(expressionData, 'VariableNames', cellstr(sampleNames));
    dataTable = addvars(dataTable, geneNames, 'Before', 1, 'NewVariableNames', 'Gene');
end

function metadata = createMetadata(sampleNames, groupA, groupB)
    groupLabels = strings(length(sampleNames), 1);

    for i = 1:length(sampleNames)
        currentName = lower(sampleNames(i));

        if contains(currentName, lower(groupA)) || contains(currentName, "patient") || contains(currentName, "disease")
            groupLabels(i) = groupA;
        elseif contains(currentName, lower(groupB)) || contains(currentName, "control") || contains(currentName, "healthy")
            groupLabels(i) = groupB;
        else
            groupLabels(i) = "Unknown";
        end
    end

    if all(groupLabels == "Unknown")
        n = length(sampleNames);
        halfPoint = floor(n / 2);
        groupLabels(1:halfPoint) = groupB;
        groupLabels(halfPoint+1:end) = groupA;
    end

    metadata = table(sampleNames(:), categorical(groupLabels), ...
        'VariableNames', {'Sample', 'Group'});
end

function X = fillMissingByGeneMean(X)
    for i = 1:size(X, 1)
        rowValues = X(i, :);
        missingValues = isnan(rowValues);

        if any(missingValues)
            rowMean = mean(rowValues, 'omitnan');
            if isnan(rowMean)
                rowMean = 0;
            end
            rowValues(missingValues) = rowMean;
            X(i, :) = rowValues;
        end
    end
end

function adjP = benjaminiHochberg(pValues)
    pValues = pValues(:);
    m = length(pValues);
    [sortedP, sortIndex] = sort(pValues, 'ascend');
    adjustedSorted = zeros(m, 1);

    for i = 1:m
        adjustedSorted(i) = sortedP(i) * m / i;
    end

    for i = m-1:-1:1
        adjustedSorted(i) = min(adjustedSorted(i), adjustedSorted(i+1));
    end

    adjustedSorted(adjustedSorted > 1) = 1;
    adjP = zeros(m, 1);
    adjP(sortIndex) = adjustedSorted;
end

function dYdt = cmlModel(~, Y, p)
    x = Y(1);
    y = Y(2);
    m = Y(3);

    dxdt = p.rL * x * (1 - x / p.KL) - p.dL * x - p.u * x;
    dydt = p.rH * y * (1 - y / p.KH) - p.alpha * x * y;
    dmdt = p.beta * x - p.gamma * m;

    dYdt = [dxdt; dydt; dmdt];
end

function [t, Y] = eulerMethod(modelFunction, timeVector, Y0)
    t = timeVector(:);
    h = t(2) - t(1);
    Y = zeros(length(t), length(Y0));
    Y(1, :) = Y0(:)';

    for i = 1:length(t)-1
        slope = modelFunction(t(i), Y(i, :)');
        Y(i+1, :) = Y(i, :) + h * slope';
        Y(i+1, :) = max(Y(i+1, :), 0);
    end
end

function [t, Y] = rk4Method(modelFunction, timeVector, Y0)
    t = timeVector(:);
    h = t(2) - t(1);
    Y = zeros(length(t), length(Y0));
    Y(1, :) = Y0(:)';

    for i = 1:length(t)-1
        currentT = t(i);
        currentY = Y(i, :)';

        k1 = modelFunction(currentT, currentY);
        k2 = modelFunction(currentT + h/2, currentY + (h/2) * k1);
        k3 = modelFunction(currentT + h/2, currentY + (h/2) * k2);
        k4 = modelFunction(currentT + h, currentY + h * k3);

        nextY = currentY + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        Y(i+1, :) = max(nextY', 0);
    end
end
