function reset_random
    % reset random numbers generator both in Matlab and Octave
    
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        rand('state',1)
        randn('state',1);
    else
        rng('default')
        %   s = RandStream('mt19937ar','Seed',2);
        %   RandStream.setGlobalStream(s);
    end
    
