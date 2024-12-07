function main
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




