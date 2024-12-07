function main
    % function main runs the robust snapshot SLAM algorithm in 45 different 
    % UE locations
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 7/12/2024
    % Tested   : Matlab version 24.1.0.2603908 (R2024a) Update 3
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    % load measurements and parameters
    [sim,params,obj] = simulation_setup;
    
    % loop through all measurement positions
    for k = 1:params.T
        dt = zeros(3,1);
        tx = sim.tx(:,k);
        y = sim.y{k};
        power = sim.power{k};

        % 1st try LoS localization
        los_candidate = find(y(1,:) == min(y(1,:)));

        [obj{k},dt(1)] = slam.localization(tx,y,power,los_candidate,params);

        if isempty(obj{k}) || sum(obj{k}.inliers) == 2 || ~los_detection(obj{k},tx,power(1,los_candidate),params)
            % NLoS condition detected --> NLoS localization
            los_candidate = [];
            [obj{k},dt(2)] = slam.localization(tx,y,power,los_candidate,params);
        end

        % mapping
        [obj{k},dt(3)] = slam.mapping(obj{k},tx,y,power,params);
        obj{k}.dt = dt;
    end
    result_summary(sim,obj,params);
end




