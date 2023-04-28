%% example showing how to convert class instances from MatLab/MTex
% for ingestion into the nomad-parser-nexus developed in the FAIRmat
% project of the German National Research Data Infrastructure
% Markus Kühbach, Humboldt-Universität zu Berlin, Department of Physics


%% initalize, specify crystal and specimen symmetries
clear;
clc;
setMTEXpref('showCoordinates','on');
setMTEXpref('FontSize',12.0);
setMTEXpref('figSize','normal');
% coordinate system, utilize SI units
% we redefine the MTex default coordinate system conventions from x2north
% and z out of plane to x east and zinto plane which is the Setting 2 case
% of TSL
% https://github.com/mtex-toolbox/mtex/issues/56
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
%right-handed Cartesian coordinate system
getMTEXpref('EulerAngleConvention');
getMTEXpref('xAxisDirection');
getMTEXpref('zAxisDirection');

%% specify test cases and their data locations
% there exists and older export_h5 function from 2013, extending this
% into an export_nexus worked smoothly until the point when one tries to
% call ipfColorKey on an ebsd class instance created inside the function
% then the [varargout{1:nargout}] = builtin('subsref',ebsd,s) cannot be
% resolved
% see publication for further details of the test datasets used

prefix = ['<<CHANGE_ME>>'];

use_case = 'ger_freiberg_hielscher';

%% currently available examples
if strcmp(use_case, 'ger_freiberg_hielscher')
    fname = [prefix '/ger_freiberg_hielscher/Forsterite.ctf'];
    CS = {... 
      'notIndexed',...
      crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
      crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite', 'color', [0.56 0.74 0.56]),...
      crystalSymmetry('12/m1', [9.7 9 5.3], [90,105.63,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside', 'color', [0.85 0.65 0.13]),...
      crystalSymmetry('m-3m', [5.4 5.4 5.4], 'mineral', 'Silicon', 'color', [0.94 0.5 0.5])};
    ebsd = EBSD.load(fname,CS,'interface','ctf', 'convertEuler2SpatialReferenceFrame');
end

disp(['Working with use_case: ' use_case]);

%% interpolate ebsd (which is an arbitrary collection of points) on square grid

disp('Interpolate EBSD data on square grid for creating NeXus default plots...');
d_x_0 = max(ebsd.unitCell(:,1)) - min(ebsd.unitCell(:,1));
d_y_0 = max(ebsd.unitCell(:,2)) - min(ebsd.unitCell(:,2));  %TODO read unit
% by construction is an axis-aligned 2D bounding box
xmin = min(ebsd.prop.x);
xmax = max(ebsd.prop.x);
ymin = min(ebsd.prop.y);
ymax = max(ebsd.prop.y);
scan_unit = ebsd.scanUnit;
if strcmp(scan_unit, 'um')
    scan_unit = lower('µm');
end

% estimate resulting size of the grid
l_x_0 = xmax - xmin;
l_y_0 = ymax - ymin;
n_x_0 = ceil(l_x_0 / d_x_0);  %TODO x and y have no meaning SUGGESTION: define where x and y is defined
n_y_0 = ceil(l_y_0 / d_y_0);

% based on this size eventually scale down the size of the image
n_max = 2048;  % one pixel to either side padding
% do not upscale smaller maps but scale larger maps to the longest side
% matching n_max
scaler = 1.;
if n_x_0 > n_max || n_y_0 > n_max
    if n_x_0 > n_y_0
        scaler = n_max / n_x_0;
    else
        scaler = n_max / n_y_0;
    end
end

% define final details of the interpolation grid
d_x = d_x_0 / scaler;
d_y = d_y_0 / scaler;
com_x = xmin + 0.5*(xmax - xmin);
com_y = ymin + 0.5*(ymax - ymin);
n_x = ceil(n_x_0 * scaler);
n_y = ceil(n_y_0 * scaler);
xmin = com_x - 0.5 * (n_x * d_x);
xmax = com_x + 0.5 * (n_x * d_x);
ymin = com_y - 0.5 * (n_y * d_y);
ymax = com_y + 0.5 * (n_y * d_y);

% interpolate ebsd using this grid
x = linspace(xmin, xmax, n_x);
y = linspace(ymin, ymax, n_y);
[x,y] = meshgrid(x,y);
xy = [x(:), y(:)].';
ebsd_interp = interp(ebsd, xy(1,:), xy(2,:));

%% write results to NeXus/HDF5
disp('Reporting results to NeXus/HDF5...');

h5w = HdfFiveSeqHdl([fname '.mtex']);
ret = h5w.nexus_create([fname '.mtex']);
ret = h5w.nexus_open('H5F_ACC_RDWR');
ret = h5w.nexus_close();

grpnm = '/entry1';
attr = io_attributes();
attr.add('NX_class', 'NXentry');
ret = h5w.nexus_write_group(grpnm, attr);

grpnm = '/entry1/indexing';
attr = io_attributes();
attr.add('NX_class', 'NXprocess');
ret = h5w.nexus_write_group(grpnm, attr);

dsnm = strcat(grpnm, '/method');
attr = io_attributes();
ret = h5w.nexus_write(dsnm, 'undefined', attr);

%% compute and add band-contrast overview image
grpnm = '/entry1/indexing/region_of_interest';
attr = io_attributes();
attr.add('NX_class', 'NXprocess');
ret = h5w.nexus_write_group(grpnm, attr);

which_descriptor = 'undefined';
if isfield(ebsd_interp.prop, 'bc')
    which_descriptor = 'normalized_band_contrast';
else
    if isfield(ebsd_interp.prop, 'confidenceindex')
        which_descriptor = 'normalized_confidence_index';
    end
end

dsnm = strcat(grpnm, '/descriptor');
attr = io_attributes();
ret = h5w.nexus_write(dsnm, which_descriptor, attr);

grpnm = '/entry1/indexing/region_of_interest/roi';
attr = io_attributes();
attr.add('NX_class', 'NXdata');
attr.add('signal', 'data');
attr.add('axes', {'axis_y', 'axis_x'});
attr.add('axis_y_indices', int64(1));
attr.add('axis_x_indices', int64(0));
ret = h5w.nexus_write_group(grpnm, attr);

dsnm = strcat(grpnm, '/data');
% compute the relevant image values ...
% the MTex-style implicit 2d arrays how they come and are used in @EBSD
if strcmp(which_descriptor, 'normalized_band_contrast')
    nxs_roi_map_u8_f = uint8(uint32( ...
        ebsd_interp.prop.bc / ...
        max(ebsd_interp.prop.bc) * 255.));
else
    nxs_roi_map_u8_f = uint8(uint32( ...
        ebsd_interp.prop.confidenceindex / ...
        max(ebsd_interp.prop.confidenceindex) * 255.));
end
% this will map NaN on zero (i.e. black in a grayscale/RGB color map)

% how it should be structured when using h5write high-level MathWorks implemented function doing internally an implicit transpose
%highlevel = uint8(zeros(fliplr([n_y n_x]))); 
%for x = 1:n_x
%    for y = 1:n_y
%        idx = y + (x-1) * n_y;
%        highlevel(x, y) = nxs_roi_map_u8_f(idx);
%    end
%end

% in contrast, how it should be structured and then flipped while writing using H5D when using the low-level MathWorks functions wrapped by HdfFiveSeqHdl
low_level = uint8(zeros([n_x n_y]));  % 731 335]));
for x = 1:n_x
    for y = 1:n_y
        idx = y + (x-1) * n_y;
        low_level(x, y) = nxs_roi_map_u8_f(idx);
    end
end

% plot matrix Matlab style as grayscale
%K = mat2gray(low_level);
%figure
%imshow(K)
%plot(ebsd_interp, ebsd_interp.bc);
%colormap gray;

attr = io_attributes();
attr.add('long_name', 'Signal');
attr.add('CLASS', 'IMAGE');
attr.add('IMAGE_VERSION', '1.2');
attr.add('SUBCLASS_VERSION', int64(15));
ret = h5w.nexus_write(dsnm, low_level, attr);

% ... and dimension scale axis positions
dsnm = strcat(grpnm, '/axis_y');
nxs_bc_y = linspace(ymin, ymax, n_y);
attr = io_attributes();
attr.add('units', scan_unit);  % TODO, convenience if larger than 1.0e or smaller than 1.e-3 auto-convert
attr.add('long_name', ['Calibrated coordinate along y-axis (', scan_unit, ')']);
ret = h5w.nexus_write(dsnm, nxs_bc_y, attr);

dsnm = strcat(grpnm, '/axis_x');
nxs_bc_x = linspace(xmin, xmax, n_x);
attr = io_attributes();
attr.add('units', scan_unit);
attr.add('long_name', ['Calibrated coordinate along x-axis (', scan_unit, ')']);
ret = h5w.nexus_write(dsnm, nxs_bc_x, attr);

dsnm = strcat(grpnm, '/title');
attr = io_attributes();
ret = h5w.nexus_write(dsnm, 'Region-of-interest overview image', attr);
% compare with what should come out
% phase_name = 'Iron fcc';
% plot(ebsd_interp(phase_name), ebsd_interp(phase_name).bc);
% ipfKey = ipfColorKey(ebsd_interp(phase_name));
% ipfKey.inversePoleFigureDirection = vector3d.Z;
% colors = ipfKey.orientation2color(ebsd_interp(phase_name).orientations);
% plot(ebsd_interp(phase_name),colors);
% colormap gray;

%% export crystal structure models
grpnm = '/entry1/indexing';
nxs_phase_id = 1;
for i = 1:length(ebsd_interp.mineralList)
    % crystallographic_database_identifier: unknown
    % crystallographic_database: unknown
    % unit_cell_abc(NX_FLOAT): 
    % unit_cell_alphabetagamma(NX_FLOAT):
    % space_group:
    % phase_identifier(NX_UINT):
    % phase_name:
    % atom_identifier:
    % atom(NX_UINT):
    % atom_positions(NX_FLOAT):
    % atom_occupancy(NX_FLOAT):
    % number_of_planes(NX_UINT):
    % plane_miller(NX_NUMBER):
    % dspacing(NX_FLOAT):
    % relative_intensity(NX_FLOAT):
    phase_name = ebsd_interp.mineralList{i};
    if strcmp(phase_name, 'notIndexed') % this is the null model
        continue;
    end
    subgrpnm = strcat(grpnm, ['/phase', num2str(nxs_phase_id)]);
    attr = io_attributes();
    attr.add('NX_class', 'NXem_ebsd_crystal_structure_model');
    ret = h5w.nexus_write_group(subgrpnm, attr);

    dsnm = strcat(subgrpnm, '/unit_cell_abc');
    unit_cell_abc = zeros([1, 3]);
    unit_cell_abc(1) = ebsd_interp.CSList{i}.aAxis.x;
    unit_cell_abc(2) = ebsd_interp.CSList{i}.bAxis.y;
    unit_cell_abc(3) = ebsd_interp.CSList{i}.cAxis.z;
    unit_cell_abc = unit_cell_abc * 0.1;  % TODO from angstroem to nm
    attr = io_attributes();
    attr.add('units', 'nm');
    ret = h5w.nexus_write(dsnm, unit_cell_abc, attr);

    dsnm = strcat(subgrpnm, '/unit_cell_alphabetagamma');
    unit_cell_alphabetagamma = zeros([1, 3]);
    unit_cell_alphabetagamma(1) = ebsd_interp.CSList{i}.alpha;
    unit_cell_alphabetagamma(2) = ebsd_interp.CSList{i}.beta;
    unit_cell_alphabetagamma(3) = ebsd_interp.CSList{i}.gamma;
    unit_cell_alphabetagamma = unit_cell_alphabetagamma / pi * 180.; % TODO from rad to deg
    attr = io_attributes();
    attr.add('units', '°');
    ret = h5w.nexus_write(dsnm, unit_cell_alphabetagamma, attr);  
    
    dsnm = strcat(subgrpnm, '/phase_identifier');
    phase_identifier = uint32(i - 1);
    attr = io_attributes();
    ret = h5w.nexus_write(dsnm, phase_identifier, attr);

    dsnm = strcat(subgrpnm, '/phase_name');
    attr = io_attributes();
    ret = h5w.nexus_write(dsnm, phase_name, attr);

    dsnm = strcat(subgrpnm, '/point_group');
    point_group = ebsd_interp.CSList{i}.pointGroup;
    attr = io_attributes();
    ret = h5w.nexus_write(dsnm, point_group, attr);

    % TODO add all the other fields relevant   

    nxs_phase_id = nxs_phase_id + 1;
end


%% export inverse pole figure (exemplified for IPF-Z) mappings and associated individual color keys
nxs_ipf_map_id = 1;
for i = 1:length(ebsd_interp.mineralList)
    phase_name = ebsd_interp.mineralList{i};
    if ~strcmp(phase_name, 'notIndexed') & sum(ebsd_interp.phaseId == i) > 0
        disp(phase_name);
        %ipf_hsv_key = ipfHSVKey(ebsd_interp(phase_name));
        ipf_key = ipfColorKey(ebsd_interp(phase_name));
        ipf_key.inversePoleFigureDirection = vector3d.Z;
        colors = ipf_key.orientation2color(ebsd_interp(phase_name).orientations);
        % from normalized colors to RGB colors
        colors = uint8(uint32(colors * 255.));
        % get the plot
        % nxs_ipf_map_u8_f = uint8(uint32(ones([n_y, n_x, 3]) * 255.));
        nxs_ipf_map_u8_f = uint8(uint32(zeros([3, n_y * n_x]) * 255.));
        nxs_ipf_y = linspace(ymin, ymax, n_y);
        nxs_ipf_x = linspace(xmin, xmax, n_x);

        % get array indices of all those pixels which were indexed as phase i
        phase_i_idx = uint32(ebsd_interp.id(ebsd_interp.phaseId == i));
        nxs_ipf_map_u8_f(:, phase_i_idx) = colors(1:length(phase_i_idx), :)';
        %for y = 1:n_y
        %    imin = 1 + (y-1)*n_x;
        %    imax = 1 + (y-1)*n_x + n_x - 1;
        %    nxs_ipf_map_dbl_f(y, :, :) = colors(tmp(imin:imax), :);
        %end
     
        grpnm = strcat('/entry1/indexing/ipf_map', num2str(nxs_ipf_map_id));
        attr = io_attributes();
        attr.add('NX_class', 'NXprocess');
        ret = h5w.nexus_write_group(grpnm, attr);

        dsnm = strcat(grpnm, '/phase_identifier');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, uint32(i - 1), attr);  
        % i-1 because in NeXus we use 0 for non-indexed but in MTex 1 and -1 in kikuchipy ...

        dsnm = strcat(grpnm, '/phase_name');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, phase_name, attr);

        dsnm = strcat(grpnm, '/projection_direction');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, single([0. 0. 1.]), attr);

        dsnm = strcat(grpnm, '/bitdepth');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, uint32(8), attr);

        dsnm = strcat(grpnm, '/program');
        attr = io_attributes();
        attr.add('version', ['Matlab: ', version ', MTex: 5.8.2']);
        ret = h5w.nexus_write(dsnm, 'mtex', attr);

%% add ipf map for specific phase
        grpnm = strcat('/entry1/indexing/ipf_map', num2str(nxs_ipf_map_id), '/ipf_rgb_map');
        attr = io_attributes();
        attr.add('NX_class', 'NXdata');
        attr.add('signal', 'data');
        attr.add('axes', {'axis_y', 'axis_x'});
        attr.add('axis_y_indices', int64(1));
        attr.add('axis_x_indices', int64(0));
        ret = h5w.nexus_write_group(grpnm, attr);

        dsnm = strcat(grpnm, '/title');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, ['Inverse pole figure color map ' phase_name], attr);

        dsnm = strcat(grpnm, '/data');
        low_level = uint8(zeros([3 n_x n_y])); %fliplr(size(nxs_ipf_map_u8_f))));
        for x = 1:n_x
            for y = 1:n_y
                idx = y + (x-1) * n_y;
                low_level(:, x, y) = nxs_ipf_map_u8_f(:, idx);
            end
        end
        attr = io_attributes();
        attr.add('long_name', 'IPF color-coded orientation mapping');
        attr.add('CLASS', 'IMAGE');
        attr.add('IMAGE_VERSION', '1.2');
        attr.add('SUBCLASS_VERSION', int64(15));
        ret = h5w.nexus_write(dsnm, low_level, attr);

        % dimension scale axis positions
        dsnm = strcat(grpnm, '/axis_y');
        % use nx_ipf_y in-place
        attr = io_attributes();
        attr.add('units', scan_unit);  % TODO, convenience if larger than 1.0e or smaller than 1.e-3 auto-convert
        attr.add('long_name', ['Calibrated coordinate along y-axis (' scan_unit ')']);
        ret = h5w.nexus_write(dsnm, nxs_ipf_y, attr);
        
        dsnm = strcat(grpnm, '/axis_x');
        % use nx_ipf_x in-place
        attr = io_attributes();
        attr.add('units', scan_unit);
        attr.add('long_name', ['Calibrated coordinate along x-axis (' scan_unit ')']);
        ret = h5w.nexus_write(dsnm, nxs_ipf_x, attr);

      
%% add specific IPF color key used
        grpnm = strcat('/entry1/indexing/ipf_map', num2str(nxs_ipf_map_id), '/ipf_rgb_color_model');
        attr = io_attributes();
        attr.add('NX_class', 'NXdata');
        attr.add('signal', 'data');
        attr.add('axes', {'axis_y', 'axis_x'});
        attr.add('axis_y_indices', int64(1));
        attr.add('axis_x_indices', int64(0));
        ret = h5w.nexus_write_group(grpnm, attr);


        dsnm = strcat(grpnm, '/title');
        attr = io_attributes();
        ret = h5w.nexus_write(dsnm, 'Inverse pole figure color key with SST', attr);

        %% get first a rendering of the color key, ...
        figure('visible','off');
        plot(ipf_key);
        % f = gcf;
        png_fnm = ['test.mtex.nxs.ipf.key.', lower(phase_name), '.png'];
        exportgraphics(gcf, png_fnm, 'Resolution', 300);
        close gcf
        % ... framegrab this image to get the pixel color values (no alpha)
        im = imread(png_fnm);
        delete(png_fnm); % remove the intermediately created figure
        % better would be to use the image directly as a matrix
        % or make a color fingerprint, i.e. defined orientation set pump
        % through color code and then rendered as n_orientations x 3 RGB
        % array this could also be used nicely for machine learning

        dsnm = strcat(grpnm, '/data');
        sz = size(im);
        low_level = uint8(zeros(fliplr(sz)));
        for x = 1:sz(2)
            for y = 1:sz(1)
                idx = y + (x-1) * n_y;
                low_level(:, x, y) = im(y, x, :);
            end
        end
        attr = io_attributes();
        attr.add('long_name', 'Signal');
        attr.add('CLASS', 'IMAGE');
        attr.add('IMAGE_VERSION', '1.2');
        attr.add('SUBCLASS_VERSION', int64(15));
        ret = h5w.nexus_write(dsnm, low_level, attr);

        % ... and dimension scale axis positions
        sz = size(im);
        dsnm = strcat(grpnm, '/axis_y');
        nxs_px_y = uint32(linspace(1, sz(1), sz(1)));
        attr = io_attributes();
        attr.add('units', scan_unit);  % TODO, convenience if larger than 1.0e or smaller than 1.e-3 auto-convert
        attr.add('long_name', 'Pixel along y-axis');
        ret = h5w.nexus_write(dsnm, nxs_px_y, attr);
        
        dsnm = strcat(grpnm, '/axis_x');
        nxs_px_x = uint32(linspace(1, sz(2), sz(2)));
        attr = io_attributes();
        attr.add('units', scan_unit);
        attr.add('long_name', 'Pixel along x-axis');
        ret = h5w.nexus_write(dsnm, nxs_px_x, attr);

        nxs_ipf_map_id = nxs_ipf_map_id + 1;
    end
end
% no path to default plottables added


%% extract MTex preferences

% plist = 'H5P_DEFAULT';
% fid = H5F.open(fnm, 'H5F_ACC_RDWR', plist); % Opens the file in read-write mode
% gid = H5G.open(fid, '/em_lab');
% gid_mtex_pref = H5G.create(gid, 'mtex_preferences', plist, plist, plist);
% H5G.close(gid);
% H5G.close(gid_mtex_pref);
% H5F.close(fid);
% 
% mtex_pref_struct = getMTEXpref;
% fn = fieldnames(mtex_pref_struct);
% for k=1:numel(fn)
%     cls = class(mtex_pref_struct.(fn{k}));
%     sz = size(mtex_pref_struct.(fn{k}));
%     disp(cls);
% %     if isstring(mtex_pref_struct.(fn{k}))
% %         disp(strcat('String\t\t'));
% %         disp(mtex_pref_struct.(fn{k}));
% %     elseif isnumeric(mtex_pref_struct.(fn{k}))
% %        disp(strcat('Numeric\t\t'));
% %        disp(mtex_pref_struct.(fn{k}));
% %     else
% %        disp('nothing');
% %     end
% % char, logical, double, function hdl, cell
% end
% subgrpnm = '/em_lab/mtex_preferences';
% h5writeatt(fnm, subgrpnm, 'xaxis_direction', mtex_pref_struct.xAxisDirection, 'TextEncoding', 'UTF-8');
% h5writeatt(fnm, subgrpnm, 'zaxis_direction', mtex_pref_struct.zAxisDirection, 'TextEncoding', 'UTF-8');
% % continue with writing other details into the file
