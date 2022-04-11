% Estimating original NA value NA_H of different datasets used in the paper
% "Minimum resolution requirements of digital pathology images for accurate 
% classification"
% (c)2019-2022 Lydia Neary-Zajiczek
% lydia.zajiczek.17@ucl.ac.uk
% Required input: dataset, Valid inputs: 'BreaKHis-4X', 'BreaKHis-10X',
% 'BreaKHis-20X', 'BreaKHis-40X', 'BACH', 'CAMELYON16-UMCU',
% 'CAMELYON16-RUMC', 'PCam'
% Optional input: cutoff_min (integer, >2, default 10) divisor of
% camera sampling frequency as lower bound of cutoff frequency range
% Optional input: show_estimation (boolean, default false)

function [] = NA_estimation(dataset, varargin)
p = inputParser;
dataset_list = {'BreaKHis4X', ...
                'BreaKHis10X', ...
            	'BreaKHis20X', ...
                'BreaKHis40X', ...
                'BACH', ...
                'CAMELYON16-UMCU', ...
                'CAMELYON16-RUMC', ...
                'PCam'};
addRequired(p,'dataset',@(x) any(strcmp(x,dataset_list)));
addOptional(p,'cutoff_min',10,@(x) isnumeric(x) && isscalar(x) && (x >= 2));
addOptional(p,'show_estimation',false,@islogical)
parse(p,dataset,varargin{:})
params = dataset_parameters(dataset);

%%%%%%%%% Location of all datasets %%%%%%%%%%%%%
base_path = 'C:\Datasets\';

% Only process benign test images for NA estimation
input_path = [base_path dataset '\data\'];
folder_list = {'test\benign\'; ...
               'test\test\'};
idx = isfolder(strcat(input_path,folder_list));
folder_list = folder_list(idx);
input_paths = strcat(input_path,folder_list);
try 
    fn_list = dir(input_paths{1});
catch
    error('Folder structure of dataset not as expected.')
end
img = imread([input_paths{1} '\' fn_list(3).name]);
img = img(:,:,1:3); %Some images had alpha channel
[ny, nx, ~] = size(img);

% Sampling parameters
d_eff = params.d_pixel/(params.M_obj*params.M_relay);
f_s_camera = 1/d_eff;
ny_scaled = params.scale_factor*ny;
nx_scaled = params.scale_factor*nx;

% Spatial frequency vectors
dx_image = d_eff/params.scale_factor;
dfx_image = 1/(dx_image*nx_scaled);
dfy_image = 1/(dx_image*ny_scaled);
fx_image = ((0:nx_scaled-1)-nx_scaled/2).*dfx_image;
fy_image = ((0:ny_scaled-1)-ny_scaled/2).*dfy_image;

% Low pass filter at f_s/2 to remove compression, aliasing artefacts
[FX_C, FY_C] = meshgrid(fx_image,fy_image);
[~, FR] = cart2pol(FX_C,FY_C);
idr = FR <= f_s_camera/2;
deg_fn_orig = double(idr);

% Range of frequencies to filter over
cutoff_freq = linspace(f_s_camera/p.Results.cutoff_min,f_s_camera/2,20);

% Pre allocate
ssim_all = zeros(length(cutoff_freq),3);
for i = 1:length(cutoff_freq)
    fn_list = dir([input_paths{1} '\*' params.img_ext]);
    if length(fn_list) > 512
        %limit size of PCam test set to 1,024 images
        fn_list = fn_list(1:256);
    end
    ssim_vec = zeros(length(fn_list),3);
    
    % Degradation function (low-pass filter, top-hat function)
    deg_fn = double(FR <= cutoff_freq(i));
    for j = 1:length(fn_list)
        img = imread([input_paths{1} '\' fn_list(j).name]);
        img_type = class(img);
        img = img(:,:,1:3); %one set of images had alpha channel
        
        % If required, upscale image
        if params.scale_factor ~= 1
            img = imresize(img,params.scale_factor,'bilinear');
        end
        if p.Results.show_estimation
            img_deg = img; 
            img_png = img;
        end
        
        % Estimate NA for each channel
        for k=1:3
            % Remove noise
            [img_in,~] = wiener2(img(:,:,k),[2 2]);
            img_ft = fft2(double(img_in));
            img_orig_ft = ifftshift(fftshift(img_ft).*deg_fn_orig);
            img_deg_ft = ifftshift(fftshift(img_ft).*deg_fn);
            switch img_type
                case 'uint8'
                    img_deg_out = uint8(abs(ifft2(img_deg_ft)));
                    img_orig_out = uint8(abs(ifft2(img_orig_ft)));
                case 'uint16'
                    img_deg_out = uint16(abs(ifft2(img_deg_ft)));
                    img_orig_out = uint16(abs(ifft2(img_orig_ft)));
            end

            % Compute similarity metric
            ssim_vec(j,k) = ssim(img_orig_out,img_deg_out);
            if p.Results.show_estimation
                img_deg(:,:,k) = img_deg_out;
                img_png(:,:,k) = img_orig_out;
            end
        end
        if p.Results.show_estimation
            if params.scale_factor ~= 1
                img_deg = imresize(img_deg,(1/params.scale_factor));
            end
            figure, imshow(img_png), title('Original Denoised Image')
            figure, imshow(img_deg), title('Low-Pass Filtered Image')
            
            % Wait for user to finish inspecting images
            input('Press enter to continue.')
        end
        disp(['Test image ' num2str(j) ' of ' num2str(length(fn_list))])
    end
    if p.Results.show_estimation
        disp(cutoff_freq(i))
        mean(ssim_vec)
    end
    disp(['Frequency ' num2str(i) ' of ' num2str(length(cutoff_freq))])
    ssim_all(i,:) = mean(ssim_vec); 
end

% Plot results to visually determine NA
figure
plot(cutoff_freq,ssim_all(:,1),'r-')
hold on
plot(cutoff_freq,ssim_all(:,2),'g-')
plot(cutoff_freq,ssim_all(:,3),'b-')
set(gca,'XDir','Reverse')
xlabel('Cutoff Frequency (1/m)')
ylabel('SSIM (a.u.)')
title(dataset)
axis([min(cutoff_freq) max(cutoff_freq) 0.7 1])

% Overlay previously determined NA_H value at green wavelength for comparison
fco = 2*params.NA_H/550e-9;
plot([fco fco],[0.7 1],'k-')
legend('Red Channel','Green Channel','Blue Channel','f_{CO} (green)')

% Save to file for dataset
saveas(gca,[base_path dataset '_NA.fig'])
save([base_path dataset '_NA.mat'])