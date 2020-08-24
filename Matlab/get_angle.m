function [crzProd, theta,costheta] = get_angle(r1,r2)

    
    norm_r1 = sqrt(r1(1)^2 + r1(2)^2);
    norm_r2 = sqrt(r2(1)^2 + r2(2)^2);
    costheta = dot(r1,r2)/(norm_r1*norm_r2);
    theta = acos(min(max(costheta,-1),1));
    crzProd = r1(1)*r2(2) - r1(2)*r2(1);
    if crzProd > 0
        theta = theta;
        crzProd = 0;
    else
        %means angle > PI (concave)
		theta = 2*pi - theta;
        crzProd = 1;
    end
    

end