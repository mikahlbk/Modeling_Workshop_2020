function [force] = get_morse_force(U, W, xsi, gamma, vector_norm,vector)

    force = (U/xsi*exp(-vector_norm/xsi)+W/gamma*exp(-vector_norm/gamma))*-vector/vector_norm;
    
end
        
  
