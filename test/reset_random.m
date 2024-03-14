function reset_random(seed)
    % reset random numbers generator both in Matlab and Octave
    
    if logical(exist('OCTAVE_VERSION', 'builtin'))
        % octave
        rand('state',seed)
        randn('state',seed);
    else
        % matlab
        rng(seed)
        %   s = RandStream('mt19937ar','Seed',seed);
        %   RandStream.setGlobalStream(s);
    end
    
