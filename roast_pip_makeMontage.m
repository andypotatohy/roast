function roast_pip_makeMontage(filename,tagname, func_mask, LeadField_tag)
% Authors: Blair Vail and Dalton Bermudez
% functions to be used after the leadField map was generated for the
% paticular subjects MRI. This pipeline is used to generate a 3x2 tDCS electrode solution
% which limits the amount of current supplied to each electrode to no more
% than 1.5 mA for a 3mA total current solution.


nii = load_untouch_nii(func_mask);% function or anatomical mask for the brain region interested in stimulation (include path to file locations and filename of the mask)
[x,y,z]= ind2sub(size(nii.img), find(nii.img==1)); %find x,y,z coordinates for RH MT
coordstim = [x,y,z];

%tagname = ['roast_targ_0087_' num2str(round) '_' num2str(cyc)];
roast_target(filename, LeadField_tag,...
    coordstim, 'coordType','voxel',...
    'optType','wls-l1per', 'orient','radial-out',...
    'desiredIntensity',0.4, 'targetRadius',10, 'elecnum', 5,... 
    'targetingTag',tagname);
ix=strfind(filename,'.');
filepath = filename(1:ix(end)-1);
% requires roast_targert modification that saves only the electrodes that apply the
% 5 maximum current in the r.montageTxt variable within the tagname.mat
% file generated (check if you have a 3X2 configuration after roast_target
% finish running)
load( [filepath '_' tagname '_targetResult.mat'])
Z = textscan(r.montageTxt,'%s','Delimiter',' ')';
for i = 1:length(Z{1})
    if strcmp(Z{1}{i}, 'mA)') || strcmp(Z{1}{i}, 'mA),')
        Z{1}{i} = '';
    elseif strcmp(Z{1}{i}(1), '(')
        Z{1}{i} =str2num(Z{1}{i}(2:end));
    end
    
    
end
emptyCells = cellfun('isempty', Z{1});
Z = {Z{1}{~emptyCells}};
orig_montage = {Z{1:2}; Z{3:4};Z{5:6};Z{7:8};Z{9:10}};
%% Scale up current to approach 0.4 V/m dosage within safety limits
minfactor = 0.4/min(r.targetMag); % scaling factor needed to bring all target voxels up to 0.4 V/m
%maxcurrent_orig = max(abs([orig_montage{:,2}])); % largest-current electrode in original montage
maxcurrent_orig = -sum(cell2mat(orig_montage(([orig_montage{:,2}]<0),2))); % total current passed through the original montage

if minfactor*maxcurrent_orig >= 3
    usefactor = 3/maxcurrent_orig; % max current injected can't exceed 4mA - if scaling up to 0.4 V/m maxes current out, use 4mA
elseif max(r.targetMag(:))*minfactor > 2.3 
    usefactor = 2.3/max(r.targetMag(:)); % field strength cannot exceed 2.3 V/m (1 scale of magnitude below most conservative estimate of tissue damage EF in Bikson et. al., 2016, Brain Stim)
else
    usefactor = minfactor; % if neither of the above conditions are met, can use scaling factor needed to bring all target voxels up to 0.4 V/m
end
% make sure the current administer by each electrode sum up to a net of
% zero (current going into the brain = current going out of the brain)
montage = orig_montage;
neg_ind = find([orig_montage{:,2}] < 0);
pos_ind = find([orig_montage{:,2}] > 0);
b = [orig_montage{:,2}];
if (sum(b(neg_ind)) == sum(b(pos_ind)))
newcurrent = usefactor*[orig_montage{:,2}]';
else
diff = abs(abs(sum(b(neg_ind))) - sum(b(pos_ind)));
min_ind = find(b == min(abs(b)));
b(min_ind) = b(min_ind) + diff;
newcurrent = usefactor*b'; 
end
neg_ind = find(newcurrent < 0);
pos_ind = find(newcurrent > 0);
diff_pos = 3 - sum(newcurrent(pos_ind)); 
diff_neg = 3 - abs(sum(newcurrent(neg_ind)));
if (sum(newcurrent(neg_ind)) == sum(newcurrent(pos_ind)))
  newcurrent = newcurrent;
elseif diff_pos ~= 0
   pos_max = find(newcurrent(pos_ind) == max(abs(newcurrent(pos_ind))));
    newcurrent(pos_max) = newcurrent(pos_max) + diff_pos;
  
elseif diff_neg ~= 0
    neg_max = find(newcurrent(neg_ind) == max(abs(newcurrent(neg_ind))));
    newcurrent(neg_max) = newcurrent(neg_max) + diff_neg;
end
% Makes sures that the current applied to each electrode does not exceed
% 1.5 mA (for patient confort) for a 3mA tDCS solution
for i = 1:length(newcurrent)
   if newcurrent(i) <= 1.5 && newcurrent(i) >= -1.5 
        newcurrent(i) = newcurrent(i);
    else
        diff_2(i) = abs(newcurrent(i)) - 1.5;
        if newcurrent(i) < 0
            newcurrent(i) = newcurrent(i) + diff_2(i);
        elseif newcurrent(i) > 0
            newcurrent(i) = newcurrent(i) - diff_2(i);
        end
    end
    
end

inx = find(diff_2 > 0);
A = newcurrent(inx);
if sum(A > 0)
pos_min = find(newcurrent == min(abs(newcurrent(pos_ind))));
newcurrent(pos_min) = newcurrent(pos_min) + diff_2(inx(A > 0));
end
if sum(A < 0)
neg_min = find(newcurrent == max(newcurrent(neg_ind)));
newcurrent(neg_min) = newcurrent(neg_min) - diff_2(inx(A < 0));
end
montage(:,2) = num2cell(newcurrent);

%% ROAST with scaled-up montage to get EF NIFTI output
recipe = {montage{1,:},montage{2,:},montage{3,:},montage{4,:},montage{5,:}};
roast(filename, recipe, 'elecType','disc','conductivities',struct('gel',4.5,'electrode',5.8e7),'resampling','on') %, 'simulationTag',...  % , 'elecSize',[10.1554 2]

end