% Two-sample t-test on smoothed, modulated, normalised GM maps (smwp1*.nii)
% HO (n = 58) vs aMCI (n = 32), nuisance covariates: TIV, age, sex.
% Regenerated from the SPM12 batch used for the analysis (Toussas et al., Brain Communications).
% Individual covariate values are read from demographics.csv (see data_templates/) rather than
% being embedded in the batch. Run with: spm('defaults','fmri'); spm_jobman('run', matlabbatch);

data_dir = fullfile(pwd, 'data');                 % <- adapt: folder containing HO/ and MCI/ subfolders
demo     = readtable(fullfile(pwd, 'demographics.csv'));   % columns: subject, group, age, sex, TIV
demo     = sortrows(demo, {'group','subject'});          % HO first, then MCI, alphabetical within group
HO  = demo(strcmp(demo.group,'HO'),:);
MCI = demo(strcmp(demo.group,'MCI'),:);

scans1 = cellfun(@(s) fullfile(data_dir,'HO', s,'mri',['smwp1' s '.nii,1']), HO.subject,  'uni',0);
scans2 = cellfun(@(s) fullfile(data_dir,'MCI',s,'mri',['smwp1' s '.nii,1']), MCI.subject, 'uni',0);

matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(pwd,'stats','2sample_TIV_age_sex_HO_vs_MCI')};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1   = scans1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2   = scans2;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept     = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca    = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova   = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c     = [HO.TIV; MCI.TIV];
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'TIV';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI  = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC   = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c     = [HO.age; MCI.age];
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI  = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC   = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c     = [HO.sex; MCI.sex];
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Sex';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI  = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC   = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {{}}, 'iCFI', {{}}, 'iCC', {{}});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Model estimation and contrasts (as run in the SPM GUI)
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(matlabbatch{1}.spm.stats.factorial_design.dir{1},'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat = matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'HO > MCI';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = 'MCI > HO';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0 0];
