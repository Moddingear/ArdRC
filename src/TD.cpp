#include "TD.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// Constructor
TD::TD(Method method, numeric_type dt, std::vector<numeric_type> params) 
    : y(0.0), yd(0.0), ydd(0.0) {
    
    if (dt <= 0) {
        throw std::invalid_argument("Sample time dt must be positive");
    }
    
    this->method = method;
    this->dt = dt;
    this->params = params;
    
    // Initialize method-specific state
    _initialize_state(0.0);
}

// Initialize state based on method
void TD::_initialize_state(numeric_type initial_value) {
    switch (method)
    {
    case Method::Euler:
        state["prev_y"] = initial_value;
        state["prev_yd"] = 0.0;
        y = initial_value;
        yd = 0.0;
        ydd = 0.0;
        break;
    case Method::Tod:
        state["v1"] = initial_value;
        state["v2"] = 0.0;
        y = initial_value;
        yd = 0.0;
        ydd = 0.0;
        break;
    case Method::LD:
        state["v1"] = initial_value;
        state["v2"] = 0.0;
        y = initial_value;
        yd = 0.0;
        ydd = 0.0;
        break;
    case Method::Red:
        state["v1"] = initial_value;
        state["v2"] = 0.0;
        y = initial_value;
        yd = 0.0;
        ydd = 0.0;
        break;
    case Method::IntD:
        state["v1"] = initial_value;
        state["v2"] = 0.0;
        y = initial_value;
        yd = 0.0;
        ydd = 0.0;
        break;
    default:
        throw std::invalid_argument("Unknown TD method: " + method);
        break;
    }
}

// Update method dispatcher
void TD::update(numeric_type ref) {
    switch (method)
    {
    case Method::Euler:
        _update_euler(ref);
        break;
    case Method::Tod:
        _update_tod(ref);
        break;
    case Method::LD:
        _update_ld(ref);
        break;
    case Method::Red:
        _update_red(ref);
        break;
    case Method::IntD:
        _update_intd(ref);
        break;
    default:
        break;
    }
}

// Reset method
void TD::reset(numeric_type initial_value) {
    _initialize_state(initial_value);
}

// Set parameters
void TD::set_parameters(std::vector<numeric_type> params) {
    this->params = params;
}

// Set sample time
void TD::set_sample_time(numeric_type dt) {
    if (dt <= 0) {
        throw std::invalid_argument("Sample time dt must be positive");
    }
    this->dt = dt;
}

// Euler-based low-pass filter with derivative estimation
void TD::_update_euler(numeric_type ref) {
    /*
     * Euler-based low-pass filter with derivative estimation.
     * Parameter: a (filter coefficient, default 0.9)
     */
    numeric_type a = !params.empty() ? params[0] : 0.9;
    
    // Low pass filter
    y = a * state["prev_y"] + (1.0 - a) * ref;
    
    // Euler approximation for first derivative
    yd = (y - state["prev_y"]) / dt;
    
    // Euler approximation for second derivative
    ydd = (yd - state["prev_yd"]) / dt;
    
    // Update state
    state["prev_y"] = y;
    state["prev_yd"] = yd;
}

// Time Optimal Differentiator
void TD::_update_tod(numeric_type ref) {
    /*
     * Time Optimal Differentiator.
     * Parameter: r (convergence rate, default 1)
     */
    numeric_type r = !params.empty() ? params[0] : 1.0;
    
    // Compute derivatives
    numeric_type dot_v1 = state["v2"];
    numeric_type dot_v2 = -r * ((state["v1"] - ref + 
                          (state["v2"] * std::abs(state["v2"])) / (2.0 * r)) > 0 ? 1.0 : -1.0);
    
    // Forward Euler integration
    state["v1"] = state["v1"] + dot_v1 * dt;
    state["v2"] = state["v2"] + dot_v2 * dt;
    
    // Outputs
    y = state["v1"];
    yd = state["v2"];
    ydd = 0.0;  // Second derivative not estimated
}

// Linear Differentiator
void TD::_update_ld(numeric_type ref) {
    /*
     * Linear Differentiator.
     * Parameter: lambda (bandwidth parameter, default 1)
     */
    numeric_type lambda_param = !params.empty() ? params[0] : 1.0;
    
    // Compute derivatives
    numeric_type dot_v1 = state["v2"] - (state["v1"] - ref) / lambda_param;
    numeric_type dot_v2 = -2.0 * (state["v1"] - ref) / lambda_param;
    
    // Forward Euler integration
    state["v1"] = state["v1"] + dot_v1 * dt;
    state["v2"] = state["v2"] + dot_v2 * dt;
    
    // Outputs
    y = state["v1"];
    yd = state["v2"];
    ydd = 0.0;  // Second derivative not estimated
}

// Robust Exact Differentiator
void TD::_update_red(numeric_type ref) {
    /*
     * Robust Exact Differentiator.
     * Parameters: lambda1, lambda2 (both default to 1)
     */
    numeric_type lambda1, lambda2;
    if (params.size() < 2) {
        lambda1 = 1.0;
        lambda2 = 1.0;
    } else {
        lambda1 = params[0];
        lambda2 = params[1];
    }
    
    // Compute derivatives
    numeric_type err = state["v1"] - ref;
    numeric_type dot_v1 = state["v2"] - lambda1 * std::sqrt(std::abs(err)) * (err > 0 ? 1.0 : -1.0);
    numeric_type dot_v2 = -lambda2 * (err > 0 ? 1.0 : -1.0);
    
    // Forward Euler integration
    state["v1"] = state["v1"] + dot_v1 * dt;
    state["v2"] = state["v2"] + dot_v2 * dt;
    
    // Outputs
    y = state["v1"];
    yd = state["v2"];
    ydd = 0.0;  // Second derivative not estimated
}

// Improved Nonlinear Tracking Differentiator
void TD::_update_intd(numeric_type ref) {
    /*
     * Improved Nonlinear Tracking Differentiator.
     * Parameters: alpha, beta, gamma, r (defaults: 1, 1, 1, 1)
     */
    numeric_type alpha, beta, gamma, r;
    if (params.size() < 4) {
        alpha = 1.0;
        beta = 1.0;
        gamma = 1.0;
        r = 1.0;
    } else {
        alpha = params[0];
        beta = params[1];
        gamma = params[2];
        r = params[3];
    }
    
    // Parameter validation
    if (alpha < 0 || alpha > 1) {
        std::cerr << "Warning: alpha should be between 0 and 1 for INTD method, clamping" << std::endl;
        alpha = std::max<numeric_type>(0.0, std::min<numeric_type>(1.0, alpha));
    }
    if (beta <= 0 || gamma <= 0 || r <= 0) {
        throw std::invalid_argument("beta, gamma, and r must be positive for INTD method");
    }
    
    // Compute derivatives
    numeric_type dot_v1 = state["v2"];
    numeric_type dot_v2 = -r * r * std::tanh((beta * state["v1"] - (1.0 - alpha) * ref) / gamma) 
                    - r * state["v2"];
    
    // Forward Euler integration
    state["v1"] = state["v1"] + dot_v1 * dt;
    state["v2"] = state["v2"] + dot_v2 * dt;
    
    // Outputs
    y = state["v1"];
    yd = state["v2"];
    ydd = 0.0;  // Second derivative not estimated
}