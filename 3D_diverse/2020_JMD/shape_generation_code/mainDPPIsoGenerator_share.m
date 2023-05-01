% ==============================================================================
% Generate only the diverse subsets of 3D isosurface unit cells.
%
% Simply run this file to generate nSamples of [res X res X res] unit cells.
% To save the unit cells, you will have to add your own code.
%
% Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020.
% "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“
% Journal of Mechanical Design. March 2021; 143(3): 031707.
%
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 2020
% ==============================================================================
dppIDs=readtable('final_dpp_sets_latent_idx_python.txt');
dppIDs = table2array(dppIDs)+1;

key = readcell('familyNameKey.csv');
% load('sfs_012_noperm_final_flat_densfilt.mat'); % only need if sampling isovalues

res = 42;
nSamples = 100;

ids = dppIDs(1,:); % row 1: joint diversity weight = 0 (shape diversity only)
n = size(ids,2);
dppGroup = zeros(n,1);
dppOrigin = zeros(n,1);
dppAB = zeros(n,1);
dppHKL = zeros(n,3);
dppLSType = cell(n,1);
dppLSField = cell(n,1);
for ii = 1:n
    % get info for structure factor generator (can also get same info from sfFunctions(i).*)
    i = ids(ii);
    kk = key{i+1,2};
    fprintf('%s\n', kk);
    ss = strsplit(kk,'-');
    dppGroup(ii) = str2double(ss{1});
    dppOrigin(ii) = str2double(ss{2});
    dppAB(ii) = str2double(ss{3});
    HKL = ss{4};
    dppHKL(ii,:) = [str2double(HKL(1)), str2double(HKL(2)), str2double(HKL(3))];
    lsType = ss{end};
    dppLSType{ii} = lsType;
    % get function handle of structure factor (same thing as sfFunctions(i).f)
    results = structureFactors(dppGroup(ii), dppHKL(ii,:), dppAB(ii), dppOrigin(ii), 1);
    f = results.f;
    % get the zero-isovalue level-set field (f<=0 is solid)
    dppLSField{ii} = getLSField2(f, 0, lsType, res, 0);
    %{
    % alternatively, generate data by sampling isovalues from feasible range
    isovals = linspace(sfFunctions(i).tRange(1), sfFunctions(i).tRange(2), nSamples+2);
    isovals = isovals(2:end-1);
    temp = cell(nSamples,1);
    % get level-set field (f<=0 is solid)
    for jj = 1:nSamples
        temp{jj} = getLSField2(f, isovals(jj), lsType, res, 0);
    end
    %}
%     figure, plotIsosurface(isosurface(dppLSField{ii}));
end


