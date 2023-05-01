function results = structureFactors(group, hkl, AB, origin, verbose)
% ==============================================================================
% Generators of level-set fields derived from crystallographic space groups,
% (currently the 36 cubic groups), normalized so all terms have coefficients = 1
% 
% f<=0 is solid
% 
% Inputs: space group number, Miller indices [h k l],
%         AB={0,1}, origin={1,2}, verbose={0,1}
% 
% Only groups 201,203,222,224,227,228 have origin=2
% Only groups 195,196,197,198,199,207,209,211,208,210,
%             212,213,214,215,216,217,218,219,220 have AB=1
%             (simplify 201-origin1, 203-origin1, 222-origin1, 224-origin1,
%              227-origin1, 227-origin1 to AB=0 only)
%
% Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020.
% "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“
% Journal of Mechanical Design. March 2021; 143(3): 031707.
% 
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 12/14/2019
% Last updated: 12/15/2019, 4/29/2020 (bug fixes)
% 
% [1] International Tables for Crystallography (2010). Vol. B, ch. 1.4, pp. 135-161.
%     https://doi.org/10.1107/97809553602060000761
% ==============================================================================

%% initialize
if ~exist('verbose', 'var'), verbose = 0; end
if ~exist('AB', 'var') || (exist('AB', 'var') && isempty(AB))
    AB = 0;
elseif AB==1 && any(~ismember(group, [195,196,197,198,199,207,209,211,208,210,...
        212,213,214,215,216,217,218,219,220])) % 4/29/20 removed 222, 224
    % if group doesn't have a B part
    results.group = group; results.f = [];
    results.exitflag = 3;
    if verbose, fprintf('    *** No. %d origin=%d *AB=%d* not allowed ***\n', group, origin, AB), end
    return
end
if ~exist('origin', 'var') || (exist('origin', 'var') && isempty(origin))
    origin = 1;
elseif origin==2 && any(~ismember(group, [201,203,222,224,227,228]))
    % if group doesn't have a second origin
    results.group = group; results.f = [];
    results.exitflag = 4;
    if verbose, fprintf('    *** No. %d *origin=%d* AB=%d not allowed ***\n', group, origin, AB), end
    return
end
if ~exist('hkl', 'var') || (exist('hkl', 'var') && isempty(hkl))
    miller = [0,1,2,3]; % allowed Miller indices components (small to limit genus of structures)
    hkl = permn(miller,3); % generate all possible permutations of the allowed components, with replacement
    hkl = hkl(2:end,:); % remove (000)
end
results = struct('group', cell(size(hkl,1),1), 'f', cell(size(hkl,1),1), ...
    'HKL', cell(size(hkl,1),1), 'AB', cell(size(hkl,1),1), ...
    'origin', cell(size(hkl,1),1), 'exitflag', cell(size(hkl,1),1));

%% get level-set functions for allowable (hkl)
c = @(u,v) cos(u*v); s = @(u,v) sin(u*v); % base functions (cos, sin)
for i = 1:size(hkl,1)
    f = [];
    h = hkl(i,1); k = hkl(i,2); l = hkl(i,3);
    results(i).group = group; results(i).AB = AB; results(i).origin = origin;
    results(i).HKL = [h k l];
    % composite functions (even, odd)
    E = @(p,q,r)@(a,b,c) ( p(h,a).*q(k,b).*r(l,c) + p(h,b).*q(k,c).*r(l,a) + p(h,c).*q(k,a).*r(l,b) );
    O = @(p,q,r)@(a,b,c) ( p(h,a).*q(k,c).*r(l,b) + p(h,c).*q(k,b).*r(l,a) + p(h,b).*q(k,a).*r(l,c) );
    switch group
        case {195,196,197} % ---------------------------------------------------
            if AB==0
                f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
            else
                f1 = E(s,s,s); f = @(x,y,z) f1(x,y,z);
            end
        case {198,199} % -----------------------------------------------------------------
            if AB==0
                if mod(h+k,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==0 && h+1==l && k==h
                    f1 = E(c,s,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l+1 && k==l
                    f1 = E(s,c,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l && k==l+1
                    f1 = E(s,s,c); f = @(x,y,z) f1(x,y,z);
                end
            else
                if mod(h+k,2)==0 && h==l && k==l
                    f1 = E(s,s,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==0 && h+1==l && k==h
                    f1 = E(s,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l+1 && k==l
                    f1 = E(c,s,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l && k==l+1
                    f1 = E(c,c,s); f = @(x,y,z) f1(x,y,z);
                end
            end
        case {200,202,204} % -----------------------------------------------------------------
            f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
        case 201 % -----------------------------------------------------------------
            if origin==1
                if mod(h+k+l,2)==0
                    f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k+l,2)==1
                    f1 = E(s,s,s); f = @(x,y,z) f1(x,y,z);
                end
            else
                if mod(h+k,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==0 && h+1==l && k==h
                    f1 = E(s,s,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l+1 && k==l
                    f1 = E(c,s,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,2)==1 && h==l && k==l+1
                    f1 = E(s,c,s); f = @(x,y,z) f1(x,y,z);
                end
            end
        case 203 % -----------------------------------------------------------------
            if origin==1
                if mod(h+k+l,4)==0
                    f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k+l,4)==1
                    f1 = E(c,c,c); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,4)==2
                    f1 = E(s,s,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k+l,4)==3
                    f1 = E(c,c,c); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            else
                if mod(h+k,4)==0 && h==l && k==l
                    f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,4)==0 && h+2==l && k==h
                    f1 = E(s,s,c); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,4)==2 && h==l+2 && k==l
                    f1 = E(c,s,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,4)==2 && h==l && k==l+2
                    f1 = E(s,c,s); f = @(x,y,z) f1(x,y,z);
                elseif mod(h+k,4)==2 && h==l && k==l
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = E(s,c,s); f4 = E(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) + f3(x,y,z) + f4(x,y,z) );
                elseif mod(h+k,4)==2 && h==l+2 && k==h
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = E(s,c,s); f4 = E(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) - f3(x,y,z) + f4(x,y,z) );
                elseif mod(h+k,4)==0 && h+2==l && k==l
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = E(s,c,s); f4 = E(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) - f3(x,y,z) - f4(x,y,z) );
                elseif mod(h+k,4)==0 && h==l && k+2==l
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = E(s,c,s); f4 = E(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) + f3(x,y,z) - f4(x,y,z) );
                end 
            end
        case {205,206} % -----------------------------------------------------------------
            if mod(h+k,2)==0 && h==l && k==l
                f1 = E(c,c,c); f = @(x,y,z) f1(x,y,z);
            elseif mod(h+k,2)==0 && h+1==l && k==h
                f1 = E(c,s,s); f = @(x,y,z) f1(x,y,z);
            elseif mod(h+k,2)==1 && h==l+1 && k==l
                f1 = E(s,c,s); f = @(x,y,z) f1(x,y,z);
            elseif mod(h+k,2)==1 && h==l && k==l+1
                f1 = E(s,s,c); f = @(x,y,z) f1(x,y,z);
            end
        case {207,209,211} % -----------------------------------------------------------------
            if AB==0
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
            else
                f1 = E(s,s,s); f2 = O(s,s,s);
                f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
            end
        case 208 % -----------------------------------------------------------------
            if AB==0
                if mod(h+k+l,2)==0
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,2)==1
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            else
                if mod(h+k+l,2)==0
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,2)==1
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            end
        case 210 % -----------------------------------------------------------------
            if AB==0
                if mod(h+k+l,4)==0
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==1
                    f1 = E(c,c,c); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,4)==2
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,4)==3
                    f1 = E(c,c,c); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            else
                if mod(h+k+l,4)==0
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,4)==1
                    f1 = E(s,s,s); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h+k+l,4)==2
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==3
                    f1 = E(s,s,s); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            end
        case 212 % -----------------------------------------------------------------
            hklmod4 = mod(h+k+l,4);
            if AB==0
                if hklmod4==0
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==1
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(c,s,s); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(s,c,s); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(s,s,c); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==2
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==3
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(c,s,s); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(s,c,s); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(s,s,c); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                end
            else
                if hklmod4==0
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==1
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(s,c,c); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(c,s,c); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(c,c,s); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==2
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==3
                    if mod(h+k,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==0 && h+1==l && k==h
                        f1 = E(s,c,c); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l+1 && k==l
                        f1 = E(c,s,c); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h+k,2)==1 && h==l && k==l+1
                        f1 = E(c,c,s); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                end
            end
        case 213 % -----------------------------------------------------------------
            hklmod4 = mod(h+k+l,4);
            if AB==0
                if hklmod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==1
                    if mod(h,2)==1 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h==k && l==h+1
                        f1 = E(c,s,s); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && k==l %h==l % BUGFIX
                        f1 = E(s,c,s); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && h==l %k==l % BUGFIX
                        f1 = E(s,s,c); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==3
                    if mod(h,2)==1 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h==k && h+1==l %h==l+1 % BUGFIX
                        f1 = E(c,s,s); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && k==l %h==l % BUGFIX
                        f1 = E(s,c,s); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && h==l %k==l % BUGFIX
                        f1 = E(s,s,c); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                end
            else % AB~=0
                if hklmod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklmod4==1
                    if mod(h,2)==1 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h==k && h+1==l %h==l+1 % BUGFIX
                        f1 = E(s,c,c); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && k==l %h==l % BUGFIX
                        f1 = E(c,s,c); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && h==l %mod(h,2)==1 && h+1==k && k==l % BUGFIX
                        f1 = E(c,c,s); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklmod4==3
                    if mod(h,2)==1 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h==k && h+1==l %h==l+1 % BUGFIX
                        f1 = E(s,c,c); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && k==l %h==l % BUGFIX
                        f1 = E(c,s,c); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && h==l %mod(h,2)==1 && h+1==k && k==l % BUGFIX
                        f1 = E(c,c,s); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                end
            end
        case 214 % -----------------------------------------------------------------
            hklMod4 = mod(h+k+l,4);
            if AB==0
                if hklMod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklMod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                end
            else
                if hklMod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklMod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                end
            end
        case {215,216,217} % -----------------------------------------------------------------
            if AB==0
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
            else
                f1 = E(s,s,s); f2 = O(s,s,s);
                f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
            end
        case {218,219} % -----------------------------------------------------------------
            hklMod2 = mod(h+k+l,2);
            if AB==0
                if hklMod2==0
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif hklMod2==1
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            else
                if hklMod2==0
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif hklMod2==1
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            end
        case 220 % -----------------------------------------------------------------
            hklMod4 = mod(h+k+l,4);
            if AB==0
                if hklMod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                elseif hklMod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(c,c,c); f2 = O(c,c,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(s,c,s); f2 = O(s,s,c);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(s,s,c); f2 = O(c,s,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(c,s,s); f2 = O(s,c,s);
                        f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                    end
                end
            else
                if hklMod4==0
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                elseif hklMod4==2
                    if mod(h,2)==0 && h==l && k==l
                        f1 = E(s,s,s); f2 = O(s,s,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==0 && h+1==k && k==l
                        f1 = E(c,s,c); f2 = O(c,c,s);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k+1 && h==l
                        f1 = E(c,c,s); f2 = O(s,c,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    elseif mod(h,2)==1 && h==k && h==l+1
                        f1 = E(s,c,c); f2 = O(c,s,c);
                        f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                    end
                end
            end
        case {221,225,229} % -----------------------------------------------------------------
            f1 = E(c,c,c); f2 = O(c,c,c);
            f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
        case 222 % -----------------------------------------------------------------
            if origin==1 % 4/29/20 moved B equations to A to simplify code
                if mod(h+k+l,2)==0 %AB==0 && 
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,2)==1 %AB~=0 && 
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            else % origin 2
                if mod(h,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==0 && h+1==k && k==l
                    f1 = E(c,s,s); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==1 && h==k+1 && h==l
                    f1 = E(s,c,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==1 && h==k && h==l+1
                    f1 = E(s,s,c); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==1 && h==k && k==l
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==1 && h==k+1 && k==l
                    f1 = E(c,s,s); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==0 && h+1==k && h==l
                    f1 = E(s,c,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==0 && h==k && h+1==l
                    f1 = E(s,s,c); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            end
        case 223 % -----------------------------------------------------------------
            if mod(h+k+l,2)==0
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
            elseif mod(h+k+l,2)==1
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
            end
        case 224 % -----------------------------------------------------------------
            if origin==1 % 4/29/20 moved B equations to A to simplify code
                if mod(h+k+l,2)==0 %AB==0 % BUGFIX: h+k+l restriction added
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,2)==1 % BUGFIX: h+k+l restriction added
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            else % origin 2
                if mod(h+k,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,2)==0 && h+1==l && k==h
                    f1 = E(s,s,c); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,2)==1 && h==l+1 && k==l
                    f1 = E(c,s,s); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,2)==1 && h==l && k==l+1
                    f1 = E(s,c,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            end
        case 226 % -----------------------------------------------------------------
            if mod(h+k+l,2)==0
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
            elseif mod(h+k+l,2)==1
                f1 = E(c,c,c); f2 = O(c,c,c);
                f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
            end
        case 227 % -----------------------------------------------------------------
            if origin==1
                if mod(h+k+l,4)==0
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==1
                    f1 = E(c,c,c); f2 = E(s,s,s); f3 = O(c,c,c); f4 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) + f3(x,y,z) - f4(x,y,z) );
                elseif mod(h+k+l,4)==2
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==3
                    f1 = E(c,c,c); f2 = E(s,s,s); f3 = O(c,c,c); f4 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) + f3(x,y,z) + f4(x,y,z) );
                end
            else % origin 2
                if mod(h+k,4)==0 && h==l && k==l %mod(h,4)==0 BUGFIX
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==0 && h+2==l && k==h %mod(h,4)==0 && h+2==k && k==l % BUGFIX
                    f1 = E(s,s,c); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l+2 && k==l %mod(h,4)==0 && h==k+2 && h==l % BUGFIX
                    f1 = E(c,s,s); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l && k==l+2 %mod(h,4)==0 && h==k && h==l+2 % BUGFIX
                    f1 = E(s,c,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l && k==l %mod(h,4)==0 && h==k && k==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) + f3(x,y,z) + f4(x,y,z) ...
                        + f5(x,y,z) + f6(x,y,z) + f7(x,y,z) + f8(x,y,z));
                elseif mod(h+k,4)==2 && h==l+2 && k==h %mod(h,4)==0 && h==k+2 && k==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) - f3(x,y,z) + f4(x,y,z) ...
                        + f5(x,y,z) - f6(x,y,z) - f7(x,y,z) + f8(x,y,z));
                elseif mod(h+k,4)==0 && h+2==l && k==l %mod(h,4)==0 && h+2==k && h==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) - f3(x,y,z) - f4(x,y,z) ...
                        + f5(x,y,z) + f6(x,y,z) - f7(x,y,z) - f8(x,y,z));
                elseif mod(h+k,4)==0 && k+2==l && h==l %mod(h,4)==0 && h==k && h+2==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) + f3(x,y,z) - f4(x,y,z) ...
                        + f5(x,y,z) - f6(x,y,z) + f7(x,y,z) - f8(x,y,z));
                end
            end
        case 228 % -----------------------------------------------------------------
            if origin==1
                if mod(h+k+l,4)==0
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==1
                    f1 = E(c,c,c); f2 = E(s,s,s); f3 = O(c,c,c); f4 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) - f3(x,y,z) - f4(x,y,z) );
                elseif mod(h+k+l,4)==2
                    f1 = E(s,s,s); f2 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k+l,4)==3
                    f1 = E(c,c,c); f2 = E(s,s,s); f3 = O(c,c,c); f4 = O(s,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) - f3(x,y,z) + f4(x,y,z) );
                end
            else % origin 2
                if mod(h+k,4)==0 && h==l && k==l %mod(h,4)==0 BUGFIX
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==0 && h+2==l && k==h %mod(h,4)==0 && h+2==k && k==l % BUGFIX
                    f1 = E(s,s,c); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l+2 && k==l %mod(h,4)==2 && h==k+2 && h==l % BUGFIX
                    f1 = E(c,s,s); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l && k==l+2 %mod(h,4)==2 && h==k && h==l+2 % BUGFIX
                    f1 = E(s,c,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h+k,4)==2 && h==l && k==l %mod(h,4)==2 && h==k && k==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) + f3(x,y,z) + f4(x,y,z) ...
                        - f5(x,y,z) - f6(x,y,z) - f7(x,y,z) - f8(x,y,z));
                elseif mod(h+k,4)==2 && h==l+2 && k==h %mod(h,4)==2 && h==k+2 && k==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) - f3(x,y,z) + f4(x,y,z) ...
                        - f5(x,y,z) + f6(x,y,z) + f7(x,y,z) - f8(x,y,z));
                elseif mod(h+k,4)==0 && h+2==l && k==l %mod(h,4)==0 && h+2==k && h==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) - f3(x,y,z) - f4(x,y,z) ...
                        - f5(x,y,z) - f6(x,y,z) + f7(x,y,z) + f8(x,y,z));
                elseif mod(h+k,4)==0 && k+2==l && h==l %mod(h,4)==0 && h==k && h+2==l % BUGFIX
                    f1 = E(c,c,c); f2 = E(c,s,s); f3 = O(s,c,s); f4 = O(s,s,c);
                    f5 = O(c,c,c); f6 = O(c,s,s); f7 = O(s,c,s); f8 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) + f3(x,y,z) - f4(x,y,z) ...
                        - f5(x,y,z) + f6(x,y,z) - f7(x,y,z) + f8(x,y,z));
                end
            end
        case 230 % -----------------------------------------------------------------
            if mod(h+k+l,4)==0
                if mod(h,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==0 && h+1==k && k==l
                    f1 = E(s,c,s); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==1 && h==k+1 && h==l
                    f1 = E(s,s,c); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                elseif mod(h,2)==1 && h==k && h==l+1
                    f1 = E(c,s,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) + f2(x,y,z) );
                end
            elseif mod(h+k+l,4)==2
                if mod(h,2)==0 && h==l && k==l
                    f1 = E(c,c,c); f2 = O(c,c,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==0 && h+1==k && k==l
                    f1 = E(s,c,s); f2 = O(s,s,c);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==1 && h==k+1 && h==l
                    f1 = E(s,s,c); f2 = O(c,s,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                elseif mod(h,2)==1 && h==k && h==l+1
                    f1 = E(c,s,s); f2 = O(s,c,s);
                    f = @(x,y,z) ( f1(x,y,z) - f2(x,y,z) );
                end
            end
        otherwise
            if verbose, fprintf('    *** invalid space group No. %d ***\n', group), end
            results(i).exitflag = 1;
            return
    end
    
    if isempty(f)
        if verbose, fprintf('    *** No. %d (%d%d%d) not allowed ***\n', group, h, k, l), end
        results(i).exitflag = 2;
    else
        results(i).f = f;
        results(i).exitflag = 0;
    end
end

end