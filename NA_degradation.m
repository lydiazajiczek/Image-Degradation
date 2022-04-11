% Degrading images from original NA value NA_H to user input value NA_L
% of different datasets used in the paper "Minimum resolution requirements 
% of digital pathology images for accurate classification"
% (c)2019-2022 Lydia Neary-Zajiczek
% lydia.zajiczek.17@ucl.ac.uk
% Required input: dataset, Valid inputs: 'BreaKHis4X', 'BreaKHis10X',
% 'BreaKHis20X', 'BreaKHis40X', 'BACH', 'CAMELYON16-UMCU',
% 'CAMELYON16-RUMC', 'PCam'
% Required input: NA_L (float, must be greater than zero)
% Optional input: save_images (boolean to save images, default true)
% Optional input: sanity_check (boolean to plot PSFs/MTFs, default false)
% Optional input: plot_ffts (boolean to plot image spectra, default false)
% Optional input: Q (integer, must be greater than zero)
% NA_L is NA value to degrade images to, Q is quantization/padding factor
% to ensure FFTs are computed correctly

function [] = NA_degradation(dataset, NA_L, varargin)
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
addRequired(p,'NA_L',@(x) isnumeric(x) && (x > 0));
addOptional(p,'save_images',true,@islogical)
addOptional(p,'sanity_check',false,@islogical)
addOptional(p,'plot_fft',false,@islogical)
addOptional(p,'Q',10,@(x) isnumeric(x) && isscalar(x) && (x >= 1));
parse(p,dataset,NA_L,varargin{:});
save_images = p.Results.save_images;
sanity_check = p.Results.sanity_check;
plot_fft = p.Results.plot_fft;
Q = p.Results.Q;

params = dataset_parameters(dataset);

%%%%%%%%% Location of all datasets %%%%%%%%%%%%%
base_path = 'C:\Datasets\';

input_path = [base_path dataset '\'];
original = 'data\';
degraded = ['data_degraded_' ...
            replace(num2str(params.NA_H,'%.2f'),'.','_') ...
            '_to_' ...
            replace(num2str(NA_L,'%.2f'),'.','_') ...
            '\'];

% train/val/test splits are different for each dataset
% Dataset may not have validation folder (PCam for example)
% Test images may have no labels (BACH for example)
folder_list = {'train\malignant\'; ...
               'train\benign\'; ...
               'val\malignant\'; ...
               'val\benign\'; ...
               'test\malignant\'; ...
               'test\benign\'; ...
               'test\test\'};

% Find which folders exist in dataset structure
idx = isfolder(strcat([input_path original],folder_list));
folder_list = folder_list(idx);
input_paths = strcat([input_path original],folder_list);
output_paths = input_paths;
try 
    fn_list = dir(input_paths{1});
catch
    error('Folder structure of dataset not as expected.')
end
img = imread([fn_list(3).folder '\' fn_list(3).name]);
img = img(:,:,1:3); %Some images had alpha channel
[ny, nx, ~] = size(img);

% Make degraded output folders here if they don't exist
if save_images
    for i = 1:length(folder_list)
        output_path = [input_path degraded folder_list{i}];
        if ~isfolder(output_path)
            mkdir(output_path)
        end
        output_paths{i} = output_path;
    end    
end

% Image processing parameters
ny_scaled = params.scale_factor*ny;
nx_scaled = params.scale_factor*nx;
lambda = [625e-9, 550e-9, 475e-9]; % From datasheet for pco edge 5.5c

% Radius of pupil functions
r_pupil_H = params.f*params.NA_H;
r_pupil_L = params.f*NA_L;

% Rayleigh criterion for minimum resolvable object separation
R_crit_H = (0.61.*lambda)./params.NA_H; 
R_crit_L = (0.61.*lambda)./NA_L;

% Incoherent cutoff frequencies
fco_H = 2*params.NA_H./lambda; 
fco_L = 2*NA_L./lambda;

% Effective pixel size in m
d_eff = params.d_pixel/(params.M_obj*params.M_relay);

% Frequency vectors of image (upsampled image size)
dx_image = d_eff/params.scale_factor;
dfx_image = 1/(dx_image*nx_scaled);
dfy_image = 1/(dx_image*ny_scaled);
fx_image = ((0:nx_scaled-1)-nx_scaled/2).*dfx_image;
fy_image = ((0:ny_scaled-1)-ny_scaled/2).*dfy_image;

% Sampling vectors in physical space (pupil plane)
Lx = Q*2*r_pupil_H; %equivalent to Q*D = n*du (field size)
Ly = (ny/nx)*Lx; %if images are not square
du = Lx/nx_scaled; 
dv = Ly/ny_scaled;
u = (du:du:Lx) - Lx/2;
v = (dv:dv:Ly) - Ly/2;
[U, V] = meshgrid(u,v);
[~,R] = cart2pol(U,V);
idx_hi = R<r_pupil_H;
idx_lo = R<r_pupil_L;

% Pupil functions
pupil_H = zeros(size(R));
pupil_H(idx_hi) = 1;
pupil_L = zeros(size(R));
pupil_L(idx_lo) = 1;

% Intensity PSFs
PSF_H = abs(fftshift(fft2(pupil_H))).^2;
PSF_L = abs(fftshift(fft2(pupil_L))).^2;
% Normalize so it sums to 1
PSF_H = PSF_H./sum(sum(PSF_H)); 
PSF_L = PSF_L./sum(sum(PSF_L)); 

% Intensity OTF/MTFs
MTF_H = abs(fftshift(fft2(PSF_H))).^2;
MTF_L = abs(fftshift(fft2(PSF_L))).^2; 
% Normalize so MTF(0,0) = 1
MTF_L = MTF_L./max(max(MTF_L));
MTF_H = MTF_H./max(max(MTF_H));

% Set to zero below machine epsilon/NaN
ep_val = 1E-9;
idneg = MTF_L < ep_val;
MTF_L(idneg) = 0;
idneg = MTF_H < ep_val;
MTF_H(idneg) = 0;

% Process all images in all folders
for i = 1:length(input_paths)
    fn_list = dir([input_paths{i} '\*' params.img_ext]);
    for j = 1:length(fn_list)
        % Read original image from dataset
        img = imread([input_paths{i} '\' fn_list(j).name]);
        img_type = class(img);
        
        % Upsample if images are very small e.g. PCam
        if params.scale_factor ~= 1
            img = imresize(img,params.scale_factor,'bilinear');
        end
        img_deg = img;
        % Fourier transform of "perfect" image convolved with PSF_H
        img_ft = fft2(double(img));
        
        % Apply channel-wise degradation function
        for k=1:length(lambda)
            % Moving from pupil to image plane is wavelength dependent
            % See Goodman 3rd ed. eqn 4-26
            % Note that images may be undersampled so MTF might be nonzero
            % over entire frequency spectrum of image!
            dx = lambda(k)/(2*params.NA_H*Q);
            x = ((0:nx_scaled-1)-nx_scaled/2).*dx; %PSF units
            dfx = (2*params.NA_H*Q)/(nx_scaled*lambda(k));
            fx = ((0:nx_scaled-1)-nx_scaled/2).*dfx; %MTF units
            dfy = (2*params.NA_H*Q)/(ny_scaled*lambda(k));
            fy = ((0:ny_scaled-1)-ny_scaled/2).*dfy; 
                        
            % Plot intensity PSF/MTF to make sure everything is sensible
            if sanity_check 
                figure %#ok<*UNRCH>
                
                % Intensity PSFs
                subplot(1,2,1)
                % High-NA (original) PSF
                plot(x,PSF_H(ny_scaled/2+1,:),'-','LineWidth',2)
                hold on
                % Low-NA (degraded) PSF
                plot(x,PSF_L(ny_scaled/2+1,:),':','LineWidth',2)
                % Rayleigh criteria
                plot([R_crit_H(j) R_crit_H(j)],[0 1],'r--','LineWidth',2) 
                plot([R_crit_L(j) R_crit_L(j)],[0 1],'g-.','LineWidth',2)
                axis([-4*R_crit_H(j) 4*R_crit_H(j) 0 max(max(PSF_H))])
                legend('PSF_H','PSF_L','d_{R,H}','d_{R,L}')
                xlabel('Image Plane Coordinates (m)')
                ylabel('Intensity PSF (a.u.)')
                
                % Intensity MTFs
                subplot(1,2,2)
                % High-NA (original) MTF
                plot(fx,MTF_H(ny_scaled/2+1,:),'-','LineWidth',2)
                hold on
                % Low-NA (degraded) MTF
                plot(fx,MTF_L(ny_scaled/2+1,:),':','LineWidth',2)
                % Incoherent cutoff frequencies
                plot([fco_H(j) fco_H(j)],[0 1],'r--','LineWidth',2)
                plot([fco_L(j) fco_L(j)],[0 1],'g-.','LineWidth',2)
                axis([-2*fco_H(j) 2*fco_H(j) 0 1])
                legend('MTF_H','MTF_L','f_{CO,H}','f_{CO,L}')
                xlabel('Image Plane Spatial Frequency (1/m)')
                ylabel('Normalized MTF (a.u.)')
                
                switch k
                    case 1
                        sgtitle('Red channel')
                    case 2
                        sgtitle('Green channel')
                    case 3
                        sgtitle('Blue channel')
                end
                
                % Wait for user to finish inspecting plot
                input('Press enter to continue.');
            end
            
            % If everything is fine, now apply to image            
            % Find overlap between fx_camera and fx_MTF
            idx = abs(fx) <= max(fx_image);
            idy = abs(fy) <= max(fy_image);
            [IDX, IDY] = meshgrid(idx,idy);
            IDZ = IDX & IDY;
            [y0,x0] = find(IDZ>0,1);
            
            % Crop MTFs to frequency spectrum of image
            MTF_H_img = imcrop(MTF_H,[x0,y0,sum(idx),sum(idy)]);
            MTF_L_img = imcrop(MTF_L,[x0,y0,sum(idx),sum(idy)]);
            MTF_H_img = imresize(MTF_H_img,size(MTF_H));
            MTF_L_img = imresize(MTF_L_img,size(MTF_L));
            deg_fn = MTF_L_img./MTF_H_img;
            idnan = isnan(deg_fn);
            deg_fn(idnan) = 0;

            % Apply to k-th channel of image FT
            img_deg_ft = fftshift(img_ft(:,:,k)).*deg_fn;
            
            % Inverse FFT
            switch img_type
                case 'uint8'
                    img_deg_out = uint8(abs(ifft2(ifftshift(img_deg_ft))));
                case 'uint16'
                    img_deg_out = uint16(abs(ifft2(ifftshift(img_deg_ft))));
            end            
            
            % Match histogram to original image to maintain color balance
            img_deg(:,:,k) = imhistmatch(img_deg_out,img(:,:,k));
            
            % See if frequency spectra look sensible before and after
            % degradation
            if plot_fft
                figure
                
                % Power spectra of images
                img_ft_p = abs(fftshift(img_ft(:,:,k)./(nx_scaled*ny_scaled)));
                img_deg_ft_p = abs(img_deg_ft(:,:,k)./(nx_scaled*ny_scaled));
                plot(fx_image,img_ft_p(ny_scaled/2,:))
                hold on
                plot(fx_image,img_deg_ft_p(ny_scaled/2,:))
                set(gca,'XLim',[min(fx_image) max(fx_image)])
                
                % Overlay MTFs and degradation function
                yyaxis right
                plot(fx_image,MTF_H_img(ny_scaled/2+1,:),'Color', ...
                     '#7E2F8E','LineWidth',2)
                hold on
                plot(fx_image,MTF_L_img(ny_scaled/2+1,:),'Color', ...
                     '#77AC30','LineWidth',2)
                plot(fx_image,deg_fn(ny_scaled/2+1,:),'LineWidth',2)
                
                % Set axis properties
                set(gca,'XLim',[min(fx_image) max(fx_image)])
                set(gca,'YLim',[0 1])
                set(gca,'YColor',[0 0 0])
                switch k
                    case 1
                        title('Red channel')
                    case 2
                        title('Green channel')
                    case 3
                        title('Blue channel')
                end
                legend('Original spectrum','Degraded spectrum','MTF_H', ...
                       'MTF_L','Degradation function')
                
                % Wait for user to finish inspecting plot
                input('Press enter to continue.');
            end
        end
            
        if params.scale_factor ~= 1
            % Downsample to original size
            img_deg = imresize(img_deg,(1/params.scale_factor));
        end
        
        % Write to disk
        if save_images
            imwrite(img_deg, [output_paths{i} '\' fn_list(j).name]);
            % Show progress
            disp(['(' num2str(j) '/' num2str(length(fn_list)) ...
                  '), folder ' num2str(i) ' of ' ...
                  num2str(length(folder_list))])
        end
    end
end
end
