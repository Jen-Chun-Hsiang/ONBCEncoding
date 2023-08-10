function angle_degrees = angle_between_vectors(vector1, vector2, type)
    if nargin == 2
        type = 'cosine';
    end
      
    % Calculate the dot product of the two vectors
    
    switch type
        case 'cosine'
            dot_product = dot(vector1, vector2);
            % Calculate the norms of the two vectors
            norm_vector1 = norm(vector1);
            norm_vector2 = norm(vector2);
            
            % Calculate the cosine of the angle between the two vectors
            cos_angle = dot_product / (norm_vector1 * norm_vector2);
            
            % Calculate the angle in radians
            angle_radians = acos(cos_angle);
            
        case 'tangent'
            dot_product = dot(vector1, vector2);
            angle_radians = atan2(det([vector1; vector2; [1 0 0]]), dot_product);
        case 'tangent2'
            angle_radians = atan2(cross(vector1,  vector2)*[1 0 0]', vector1*vector2');
    end
    % Convert the angle to degrees
    angle_degrees = rad2deg(angle_radians);
end