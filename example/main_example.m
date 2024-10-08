% ---------------------------------------------------------------------
%   HOUSEKEEPING
% ---------------------------------------------------------------------
clc;
clear;
close all;
warning('off', 'all')   

%Set seed
rng(14);

addpath(genpath('../src')); % This contains the functions for bootstrap correction
addpath(genpath('./CMG')); % This contains the functions for building the preconditioner for Laplacian systems


% ---------------------------------------------------------------------
%   SOURCES AND PARAMETERS
% ---------------------------------------------------------------------   

% Bootstrap parameters    
N_boot = 100; 
n_boot = N_boot;
n_lev = N_boot;
mytol = 1e-6;
ind_export = 1;
ind_light = 0;
LdM_mom = 0;

% Source data
nm_data = '.\data\example_large.csv';    

% ----------------------------------------------------------------------
%   LOAD DATA
% ----------------------------------------------------------------------

%% READ DATA
data=readtable(nm_data);

id=data.id;
firmid=data.firmid;
pe_t_0=data.pe_t;
fe_t_0=data.fe_t;
age = data.age;
age_sq = age.^2;
y = data.y;
year = data.year;   
clear data

% Create a cluster identifier (match)
cluster = findgroups(id,firmid); %cluster at the match level

% Fixed effects, controls
mat_controls = age_sq;
other_fe = year;
clear age age_sq

% Create different grouping vectors
group1 = randi([1 20],1,size(y,1))';
group2 = randi([1 50],1,size(y,1))';
group3 = randi([1 800],1,size(y,1))';

group_mat = [group1 group2 group3];
v_filename_group = {'1'; '2' ; '3'};


% Start paralell environment
try
    pool=parpool;
end
%% Tests    
disp('------------------------------------')
disp('  1. MINIMAL EXAMPLE WITHOUT PRINTING')
disp('------------------------------------')

rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls);

time_minimal = toc(start);



disp('------------------------------------')
disp('  2. MINIMAL EXAMPLE PRINTING')
disp('------------------------------------')

filename = strcat('.\results\minimal_overall');

rng(435)
start = tic;
correction_all(y, id, firmid, other_fe, mat_controls, [], [], [], ...
    [], [], [], [], [], [], [],...
    filename, ind_export, []);

time_minimal_print = toc(start);


disp('------------------------------------')
disp('  3. EXAMPLE HCU')
disp('------------------------------------')

type_hc = 'hc_u';
type_leave = 'match';
filename = strcat('.\results\',type_hc,'_',type_leave,'_overall');



rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, [], [], [], ...
    type_hc, [], [], [], [], [], [],...
    filename, ind_export, []);

time_hcu = toc(start);


disp('------------------------------------')
disp('  4. EXAMPLE HCU MULTI GROUP')
disp('------------------------------------')

type_hc = 'hc_u';
type_leave = 'match';
group = group_mat;
filename = strcat('.\results\',type_hc,'_',type_leave,'_per_group');


rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, group, [], [], ...
    type_hc, [], [], [], [], [], [],...
    filename, ind_export, []);

time_hcu_multi = toc(start);


disp('------------------------------------')
disp('  5. EXAMPLE CLUSTER-OUT')
disp('------------------------------------')

type_hc = 'hc_u_clus';
type_leave = 'cluster';
filename = strcat('.\results\',type_hc,'_',type_leave,'_overall');

rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, [], [], [], ...
    type_hc, [], cluster, [], [], [], [],...
    filename, ind_export, []);

time_cluster = toc(start);

disp('------------------------------------')
disp('  6. EXAMPLE CLUSTER-OUT MULTI GROUP')
disp('------------------------------------')

type_hc = 'hc_u_clus';
group = group_mat;
type_leave = 'cluster';
filename = strcat('.\results\',type_hc,'_',type_leave,'_per_group');


rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, group, [], [], ...
    type_hc, [], cluster, [], [], [], [],...
    filename, ind_export, []);


time_cluster_multi = toc(start);


disp('------------------------------------')
disp('  7. EXAMPLE MATCH-OUT WITH LdM MOMENT')
disp('------------------------------------')

type_hc = 'hc_u_match';
group = [];
type_leave = 'match';
filename = strcat('.\results\',type_hc,'_',type_leave,'_overall');


rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, [], [], [], ...
    [], [], [], [], [], LdM_mom, year,...
    filename, ind_export, []);

time_match = toc(start);



disp('------------------------------------')
disp('  9. CHECKING SAMPLE SELECTION')
disp('------------------------------------')

rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, [], ...
    [], [], 'hc_u', 'match',...
    [], ind_light, [], [], [],...
    [], [], []);

time_sample_sel = toc(start);



disp('------------------------------------')
disp('  10. FULLY FLEXIBLE')
disp('------------------------------------')

filename = strcat('.\results\flexible');

rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, group,...
			n_lev, n_boot, type_hc, type_leave,...
			cluster, ind_light, mytol, LdM_mom, year,...
			filename, ind_export, v_filename_group);

time_flexible = toc(start);


disp('------------------------------------')
disp('  11. MEMORY EFFICIENT WITHOUT PRINTING')
disp('------------------------------------')

filename = strcat('.\results\memory');

rng(435)
start = tic;
[plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, [], [], [], ...
    [], [], [], ind_light, [], [], [],...
    [], [], []);

time_memory = toc(start);

