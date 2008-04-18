function plotipdf(odf,r,varargin)
% plot inverse pole figures
%
%% Input
%  odf - @ODF
%  r   - @vector3d specimen directions
%
%% Options
%  RESOLUTION - resolution of the plots
%  
%% Flags
%  REDUCED  - reduced pdf
%  COMPLETE - plot entire (hemi)-sphere
%
%% See also
% S2Grid/plot savefigure

argin_check(r,'vector3d');

% plotting grid
h = S2Grid('PLOT',varargin{:},...
  'MAXTHETA',rotangle_max_y(odf(1).CS,varargin{:})/2,...
  'MAXRHO',max(2*pi*check_option(varargin,'COMPLETE'),...
  rotangle_max_z(odf(1).CS,varargin{:})));

multiplot(@(i) h,@(i) pdf(odf,h,r(i)./norm(r(i)),varargin{:}),length(r),...
  'DISP',@(i,Z) [' iPDF r=',char(r(i)),' Max: ',num2str(max(Z(:))),...
  ' Min: ',num2str(min(Z(:)))],...
  'ANOTATION',@(i) char(r(i),'LATEX'),...
  'MINMAX','SMOOTH',...
  varargin{:});

set(gcf,'Name',['recalculated Pole Figures of Specimen ',inputname(1)]);
