function flag = los_detection(obj,tx,power,params)
    % perform statistical test based on the negative log-likeihood to check
    % is the candidate LoS path in line with the log distance path-loss
    % model
    
    if isempty(obj)
        flag = false;
        return
    else
        theta =  params.pathloss.theta;
        threshold = params.pathloss.threshold;
        
        d = sqrt(sum((obj.xn(1:2) - tx(1:2)).^2));
        nu = (theta(1) + theta(2)*10*log10(d)) - power;
        
        L = -0.5*(log(2*pi) + log(theta(3)) + nu^2/theta(3));
        if  L < threshold
            flag = false;
            return
        end
    end
    flag = true;
end