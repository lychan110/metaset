% ==============================================================================
% Generate the level-set fields and voxels of 3D isosurface unit cells.
%
% Simply run this file to generate nSamples of [res X res X res] unit cells.
% To save the unit cells, you will have to add your own code.
% 
% Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020. 
% "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“ 
% Journal of Mechanical Design. March 2021; 143(3): 031707.
%
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 5/2/2019
% Late updated: 5/11/2019, 6/18/2019
% ==============================================================================
key = readcell('familyNameKey.csv');
load('sfs_012_noperm_final_flat_densfilt.mat');

res = 64;
nSamples = 100;
n = length(sfFunctions);
skipid = [112,130,157,189,200,211,236,270,283,302];

for ii = 1:n
    if ~ismember(ii,skipid)
        % get info for structure factor generator (can also get same info from sfFunctions(i).*)
        kk = key{ii+1,2};
        fprintf('%s\n', kk);
        ss = strsplit(kk,'-');
        spaceGroup = str2double(ss{1});
        origin = str2double(ss{2});
        AB = str2double(ss{3});
        HKL = ss{4};
        HKL = [str2double(HKL(1)), str2double(HKL(2)), str2double(HKL(3))];
        lsType = ss{end};
        % get function handle of structure factor (same thing as sfFunctions(i).f)
        results = structureFactors(spaceGroup, HKL, AB, origin, 1);
        f = results.f;
        % get the zero-isovalue level-set field (f<=0 is solid)
        levelSetField = getLSField2(f, 0, lsType, res, 0);
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
%         figure, plotIsosurface(isosurface(levelSetField));
        %{
        % get density (volume fraction)
        vox = (levelSetField<=0);
        dens = mean(vox(:));
%         figure, plotIsosurface(vox);
        %}
    end
end


