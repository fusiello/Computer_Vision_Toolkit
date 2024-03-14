
%restoredefaultpath;

root  = pwd;
addpath(fullfile(root,'m-files'             ));
addpath(fullfile(root,'m-files','aux_fun'       ));
addpath(fullfile(root,'thirdparty'          ),'-end');

savepath;

if ~logical(exist('OCTAVE_VERSION', 'builtin'))
    help m-files
else
    what('m-files')
end

