#pragma once

#include <vector>
#include <string>
#include <map>

class TD {
public:
    typedef float numeric_type;
    enum class Method {
        Euler, //Euler-based low passs filter with derivative estimation
        Tod, //Time Optimal Differentiator
        LD, //Linear Differentiator
        Red, //Robust Exact Differentiator
        IntD, //Improved nonlinear tracking differentiator
    };
    // Output properties
    numeric_type y;      // Filtered output
    numeric_type yd;     // First derivative
    numeric_type ydd;    // Second derivative
    
    // Constructor
    TD(Method method, numeric_type dt, std::vector<numeric_type> params = {});
    
    // Update method
    void update(numeric_type ref);
    
    // Reset method
    void reset(numeric_type initial_value = 0.0);
    
    // Set parameters
    void set_parameters(std::vector<numeric_type> params);
    
    // Set sample time
    void set_sample_time(numeric_type dt);

private:
    Method method;
    numeric_type dt;
    std::vector<numeric_type> params;
    std::map<std::string, numeric_type> state;
    
    // Initialization
    void _initialize_state(numeric_type initial_value);
    
    // Update methods for each differentiator type
    void _update_euler(numeric_type ref);
    void _update_tod(numeric_type ref);
    void _update_ld(numeric_type ref);
    void _update_red(numeric_type ref);
    void _update_intd(numeric_type ref);
};