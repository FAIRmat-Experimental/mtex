function display(SO3F,varargin)
% standard output

if check_option(varargin,'skipHeader')
  disp('  <strong>function handle component</strong>');
else
  displayClass(SO3F,inputname(1),[],'moreInfo',symChar(SO3F),varargin{:});
  if numel(SO3F) > 1, disp(['  size: ' size2str(SO3F)]); end
end

try
  m = mean(SO3F,'all','bandwidth',16);
  disp(['  weight: ' xnum2str(m(1))])
end
disp(' ')

end