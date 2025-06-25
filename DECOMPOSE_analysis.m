clear; clc; close all;
%% THIS IS MY INTERPRETATION OF DECOMPOSE QSM
% Code by Tim Ho (timothy.ho@virginia.edu)


% MEDI toolbox (http://pre.weill.cornell.edu/mri/pages/qsm.html)
%addpath(genpath("D:\SOFTWARE\QSM\MEDI_toolbox"))

% Set QSM tool directory path 
% STI Suite (Version 3.0) (https://people.eecs.berkeley.edu/~chunlei.liu/software.html)
%addpath(genpath("D:\SOFTWARE\QSM\STISuite_V3.0"))

% Set x-separation tool directory path
% CHI-SEP Toolbox (https://github.com/SNU-LIST/chi-separation)
%addpath(genpath("D:\SOFTWARE\Chisep_Toolbox_v1.2"))

addpath(genpath("utils"))

%% HEADER VARIABLES
CF = 123200000;
B0 = 2.89;               % T
B0_direction = [0 0 1];
TE = (1:8) * 3.5 * 1e-3; % Echo times in s
TE = TE(1:end);
delta_TE = TE(2)-TE(1);
voxel_size = [1 1 1];
gamma = 42.58 * 2*pi; % rad/T

%% Open data sets
folder = 'D:\Recon\ANALYSIS\PHANTOM_FE3O4_CACO3_06172025';
mag_multi_echo = niftiread([folder '\' 'img.nii']);
phs_multi_echo = niftiread([folder '\' 'phantom\results\combined_phase.nii']);
mask_brain = double(niftiread([folder '\Segmentation.nii'])); % [X Y Z]
mag_multi_echo = mag_multi_echo/max(mag_multi_echo(:)); % Normalize
mag_multi_echo = squeeze(sqrt(sum(mag_multi_echo.^2,5)));
mag_multi_echo = double(mag_multi_echo);
phs_multi_echo = double(squeeze(phs_multi_echo));
mag_multi_echo = mag_multi_echo(:,:,:,:);
phs_multi_echo = phs_multi_echo(:,:,:,:);

%EXPECTED INPUT for both Magnitude image and phase image [X Y Z TE]
[N1,N2,N3,N4]  = size(mag_multi_echo);

%% QSM
% Normal QSM
[phase, N_std] = Preprocessing4Phase(mag_multi_echo, phs_multi_echo);
pad_size=[12 12 12];
[field_map, ~] = MRPhaseUnwrap(phase,'voxelsize',voxel_size,'padsize',pad_size);
smv_size=25;
[local_field, mask_brain_new]=V_SHARP(field_map, mask_brain,'voxelsize',voxel_size,'smvsize',smv_size);
pad_size = [12, 12, 12];
QSM = QSM_star(local_field,mask_brain_new,'TE',delta_TE*1e3,'B0',B0,'H',B0_direction,'padsize',pad_size,'voxelsize',voxel_size);
niftiwrite(QSM,[folder '\' 'QSM.nii.gz'])

%% Individual QSM for Each Echo
QSM_multi = zeros(N1,N2,N3);
for i = 1:length(TE)
    % Phase Unwrapping
    % Laplacian-based method from STI Suite
    pad_size=[12 12 12];
    phase = squeeze(phs_multi_echo(:,:,:,i));
    [field_map, ~] = MRPhaseUnwrap(phase,'voxelsize',voxel_size,'padsize',pad_size);
    field_map_hz = field_map / (2*pi*TE(i)); % convert rad to hz
    
    % Background field removal
    % V-SHARP from STI Suite
    smv_size=25;
    [local_field, mask_brain_new]=V_SHARP(field_map, mask_brain,'voxelsize',voxel_size,'smvsize',smv_size);
    local_field_hz = local_field / (2*pi*TE(i)); % convert rad to hz
    
    % QSM
    % 1. STARQSM from STI Suite
    pad_size = [12, 12, 12];
    QSM_multi(:,:,:,i) = QSM_star(local_field,mask_brain,'TE',TE(i)*1e3,'B0',B0,'H',B0_direction,'padsize',pad_size,'voxelsize',voxel_size);
end
save('QSM',"QSM","QSM_multi","mask_brain","TE","B0")


%% DECOMPOSE DATA
y_data = zeros(N1,N2,N3,length(TE));
for i = 1:length(TE)
    mag = mag_multi_echo(:,:,:,i); %S0_map.*exp(-TE(i)/(T2_map*1e-3));
    y_data(:,:,:,i) = (mag.*exp(-1i*(2/3).*squeeze(QSM_multi(:,:,:,i))*gamma*B0*TE(i))).*mask_brain_new;
end
save('y_data',"y_data", "N1", "N2", "N3", "TE")

%% Fit Model
% Change varnames later
C_plus_map    = zeros(N1,N2,N3);
C_minus_map   = zeros(N1,N2,N3);
C0_map        = zeros(N1,N2,N3);
R0_map        = zeros(N1,N2,N3);
chi_plus_map  = zeros(N1,N2,N3);
chi_minus_map = zeros(N1,N2,N3);

%% OPTIMIZATION PARAMETERS
N_inner       = 10;

% Set optimization options for lsqcurvefit
% 'UseParallel' set to false since data if fitting per voxel
% 'Display', 'iter'
%options1 = optimoptions('lsqcurvefit', 'Display', 'off', 'UseParallel', false, 'SpecifyObjectiveGradient', true);
options1 = optimoptions('lsqcurvefit', 'Display', 'off', 'UseParallel', false);
options2 = optimoptions('lsqcurvefit', 'Display', 'off', 'UseParallel', false);
options3 = optimoptions('lsqcurvefit', 'Display', 'off', 'UseParallel', false);
% Define the bounds for the variables
lb  = [0, 0, 0];
ub  = [inf, inf, inf]; %[1, 1, 1]; % BK Bounds %[inf, inf, inf]; % JH Bounds for 5.1
lb2 = [0];
ub2 = [inf];
upper = 0.5;
lb3 = [0, -upper];
ub3 = [upper,  0];       % JH Bounds for 5.3 [0.1,0.1] 11.4T; Berkley was [0.5,0.5] 3T
tic
parfor i =1:N1 % Parallelize Slice
    fprintf('Slice: %.3d \n', i);
    for j = 1:N2
        for k = 1:N3 % 1:N3 Everything
            ydata = squeeze(y_data(i,j,k,:))';
            if mask_brain(i,j,k) == 0 || sum(ydata) == 0
                continue
            end
            %fprintf('Slice: %.3d %.3d %.3d\n', i, j, k)
           
            % Initial guess for the variables chi_plus, chi_minus
            chi_0 = [0.05, -0.05]; % ppm 
            % Initial guess for R_star_2_0
            R0 = 25;              % Hz         
            % Initial guess for the variables C_plus, C_minus, C_0
            C0 = [0.3, 0.3, 0.4]; % Partial Fractions

            for iteration = 1:N_inner
                %% Optimization 1
                % Fit Variables Cplus, Cminus, C0
                % Fixed Variables Chi_plus,Chi_minus,R0
                chi_plus  = chi_0(1);
                chi_minus = chi_0(2);
                % Solve the non-linear least squares problem using lsqcurvefit
                model1 = @(x,x_data) complex_to_real(objective1(x, x_data, chi_plus, chi_minus, R0, B0));
                [C0, resnorm, residual, exitflag, output] = lsqcurvefit(model1, C0, TE, complex_to_real(ydata), lb, ub, options1);
                
                %% Optimization 2
                % Fit Variables R0
                % Replace Cplus, Cminus, C0 variables
                % Fixed Variables Chi_plus,Chi_minus
                C_plus  = C0(1);
                C_minus = C0(2);
                C_0     = C0(3);
                model2 = @(x,x_data) complex_to_real(objective2(x, x_data, C_plus, C_minus, C_0, chi_plus, chi_minus, B0));
                [R0, resnorm, residual, exitflag, output] = lsqcurvefit(model2, R0, TE, complex_to_real(log(ydata)), lb2, ub2, options2);

                %% Optimization 3
                % Fit chi_plus, chi_minus
                % Replace R0 (completed in Optimization 2)
                % Fixed Variables Cplus, Cminus, C0
                model3 = @(x,x_data) complex_to_real(objective3(x, x_data, C_plus, C_minus, C_0, R0, B0));
                [chi_0, resnorm, residual, exitflag, output] = lsqcurvefit(model3, chi_0, TE, complex_to_real(log(ydata)), lb3, ub3, options3);

            end

            %% Save Data
            C_plus_map(i,j,k)    = C0(1);
            C_minus_map(i,j,k)   = C0(2);
            C0_map(i,j,k)        = C0(3);
            R0_map(i,j,k)        = R0;
            chi_plus_map(i,j,k)  = chi_0(1);
            chi_minus_map(i,j,k) = chi_0(2);

        end
    end
end
toc
save('DECOMPOSE', "C_plus_map","C_minus_map", "C0_map", "chi_plus_map","chi_minus_map","R0_map")

%% Create Composites
gamma = 42.58 * 2*pi;
a = (2*pi*gamma*B0)/(9*sqrt(3));
den = (2/3) * gamma * B0 * sum(TE);
PSC = zeros(size(mask_brain));
DSC = zeros(size(mask_brain));
COM = zeros(size(mask_brain));

% Compute susceptibility for each individual pixel
for i = 1:N1
    for j = 1:N2
        for k = 1:N3
            if mask_brain(i,j,k) == 0
                continue
            end
chi_pos = pscModel(TE, C_plus_map(i,j,k), C_minus_map(i,j,k), C0_map(i,j,k), chi_plus_map(i,j,k), chi_minus_map(i,j,k), R0_map(i,j,k), B0);
chi_pos = sum(angle(chi_pos))/den;
chi_neg = dscModel(TE, C_plus_map(i,j,k), C_minus_map(i,j,k), C0_map(i,j,k), chi_plus_map(i,j,k), chi_minus_map(i,j,k), R0_map(i,j,k), B0);
chi_neg = sum(angle(chi_neg))/den;
chi_tot = signalModel(TE, C_plus_map(i,j,k), C_minus_map(i,j,k), C0_map(i,j,k), chi_plus_map(i,j,k), chi_minus_map(i,j,k), R0_map(i,j,k), B0);
chi_tot = sum(angle(chi_tot))/den;
PSC(i,j,k) = -chi_pos;
DSC(i,j,k) = -chi_neg;
COM(i,j,k) = chi_tot;
        %disp(['NEG ', num2str(chi_neg)])
        %disp(['POS ', num2str(-chi_pos)])
        end
    end
end

niftiwrite(DSC,'XPSC')
niftiwrite(PSC,'XDSC')
niftiwrite(COM,'XCOM')

