function ang = meanAngle(inputangles)
    inputangles = inputangles(:);
    inputangles = inputangles(~isnan(inputangles));
    x = mean(cos(inputangles*pi/180));
    y = mean(sin(inputangles*pi/180));
%     keyboard;
    if abs(x) < 1e-6
        ang = 90;
    else
        ang = atan2(y, x)*180/pi;
    end
end