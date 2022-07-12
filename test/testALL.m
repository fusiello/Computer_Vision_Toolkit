close all
clear functions

reset_random;

CVT = fileparts(pwd);

testTWO
testREG
testMULTI
testCARVE
testIMG
testAUX

% find called mfiles
mfiles = inmem('-completenames');
called ={};

% select the ones in home dir
for i =1:size(mfiles,1)
    if  ~isempty(findstr(mfiles{i}, CVT))
        %disp(mfiles{i});
        called = [called ; mfiles{i}];
    end
end

all_files  = getAllFiles(CVT,{'m'});

disp('These m-files are not covered by test:')
foo=setdiff(all_files, called);
fprintf('%s\n',foo{:});