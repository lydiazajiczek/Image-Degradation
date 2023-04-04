% Function that returns necessary imaging system parameters for each of the
% datasets used in the paper "Minimum resolution requirements of digital 
% pathology images for accurate classification"
% (c)2019-2022 Lydia Neary-Zajiczek
% lydia.zajiczek.17@ucl.ac.uk
%
% Dataset-specific parameters. Valid inputs: 'BreaKHis4X', 'BreaKHis10X',
% 'BreaKHis20X', 'BreaKHis40X', 'BACH', 'CAMELYON16_UMCU',
% 'CAMELYON16_RUMC', 'PCam'
% Returns 'params' struct containing:
% 'NA_H': Original NA of imaging system
% 'M_obj': Magnification of objective
% 'M_relay': Magnification of any relay lens present, 1 otherwise 
% 'dx_pixel': Physical size of camera pixel in m 
% 'f': Focal length of tube lens in m  
% 'img_ext': Extension of image files, e.g. '.png', '.tif'
            
function params = dataset_parameters(dataset)

% Images from PCam must be upsampled for frequency-space manipulations to
% work properly, set to 1 for all other datasets
scale_factor = 1;
switch dataset
    case 'BACH' %Leica DM 2000 LED microscope with Leica ICC50 HD Camera
        NA_H = 0.30; 
        M_obj = 20;  
        M_relay = 0.5; 
        d_pixel = 3.2E-6; 
        f = 200e-3/M_obj; %Leica, Nikon, Mitutoyo tube lens fl 200 mm  
        img_ext = '.tif';
    case 'CAMELYON16-UMCU' %NanoZoomer-XR, for NA estimation only
        NA_H = 0.75; 
        M_obj = 20;
        M_relay = 2;
        d_pixel = 9.2E-6;
        f = 200e-3/M_obj; %Leica, Nikon, Mitutoyo tube lens fl length 200 mm
        img_ext = '.tiff';
    case 'CAMELYON16-RUMC' %Pannoramic 250 Flash II, for NA estimation only
        NA_H = 0.80; %From datasheet
        M_obj = 20;  
        M_relay = 1; 
        d_pixel = 5.5E-6;
        f = 200e-3/M_obj; %Leica, Nikon, Mitutoyo tube lens fl 200 mm 
        img_ext = '.tiff';
    case 'BreaKHis4X' %Olympis BX-50 w extended apochromats
        NA_H = 0.16;
        M_obj = 4;
        M_relay = 3.3;
        d_pixel = 6.5E-6;
        f = 180e-3/M_obj; %Olympus tube lens fl 180 mm
        img_ext = '.png';
    case 'BreaKHis10X'
        NA_H = 0.40;
        M_obj = 10;
        M_relay = 3.3;
        d_pixel = 6.5E-6;
        f = 180e-3/M_obj; %Olympus tube lens fl 180 mm
        img_ext = '.png';
    case 'BreaKHis20X'
        NA_H = 0.80;
        M_obj = 20;
        M_relay = 3.3;
        d_pixel = 6.5E-6;
        f = 180e-3/M_obj; %Olympus tube lens fl 180 mm
        img_ext = '.png';
    case 'BreaKHis40X'
        NA_H = 1.4;
        M_obj = 40;
        M_relay = 3.3;
        d_pixel = 6.5E-6;
        f = 180e-3/M_obj; %Olympus tube lens fl 180 mm
        img_ext = '.png';
    case 'PCam'
        NA_H = 0.13;
        M_obj = 10;
        M_relay = 1;
        d_pixel = 9.72E-6; 
        f = 200e-3/M_obj; %Leica, Nikon, Mitutoyo tube lens fl 200 mm 
        img_ext = '.tif';
        scale_factor = 10;
    case 'BRACS'
        NA_H = 0.75;
        M_obj = 20;
        M_relay = 2;
        d_pixel = 10E-6;
        f = 200e-3/M_obj; %Leica, Nikon, Mitutoyo tube lens fl 200 mm 
        img_ext = 'png';
    otherwise
        error('Unknown dataset. See function help for valid datasets.')
end

params = struct('NA_H',NA_H,'M_obj',M_obj,'M_relay',M_relay,  ...
                'd_pixel',d_pixel,'f',f,'img_ext',img_ext, ...
                'scale_factor',scale_factor);
end