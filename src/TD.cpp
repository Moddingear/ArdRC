#include "TD.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// Constructor
TD::TD(numeric_type dt) 
{
	
	if (dt <= 0) {
		throw std::invalid_argument("Sample time dt must be positive");
	}

	this->dt = dt;
	
	// Initialize state
	reset(0.0);
}

// Reset method
void TD::reset(numeric_type initial_value) {
	y = initial_value;
	yd = 0.0;
	ydd = 0.0;
}

// Set sample time
void TD::set_sample_time(numeric_type dt) {
	if (dt <= 0) {
		throw std::invalid_argument("Sample time dt must be positive");
	}
	this->dt = dt;
}

TD_Euler::TD_Euler(TD::numeric_type dt, TD::numeric_type a)
	: TD(dt)
{
	this->a = a;
}

// Euler-based low-pass filter with derivative estimation
void TD_Euler::update(numeric_type ref) {
	/*
	 * Euler-based low-pass filter with derivative estimation.
	 * Parameter: a (filter coefficient, default 0.9)
	 */
	
	// Low pass filter
	y = a * prev_y + (1.0 - a) * ref;
	
	// Euler approximation for first derivative
	yd = (y - prev_y) / dt;
	
	// Euler approximation for second derivative
	ydd = (yd - prev_yd) / dt;
	
	// Update state
	prev_y = y;
	prev_yd = yd;
}

void TD_Euler::reset(numeric_type initial_value) {
	prev_y = initial_value;
	prev_yd = 0.0;
	TD::reset(initial_value);
}

TD_TOD::TD_TOD(numeric_type dt, numeric_type r)
	: TD(dt)
{
	this->r = r;
}

// Time Optimal Differentiator
void TD_TOD::update(numeric_type ref) {
	/*
	 * Time Optimal Differentiator.
	 * Parameter: r (convergence rate, default 1)
	 */
	
	// Compute derivatives
	numeric_type dot_v1 = v2;
	numeric_type dot_v2 = -r * ((v1 - ref + 
						  (v2 * std::abs(v2)) / (2.0 * r)) > 0 ? 1.0 : -1.0);
	
	// Forward Euler integration
	v1 = v1 + dot_v1 * dt;
	v2 = v2 + dot_v2 * dt;
	
	// Outputs
	y = v1;
	yd = v2;
	ydd = 0.0;  // Second derivative not estimated
}

void TD_TOD::reset(numeric_type initial_value) {
	v1 = initial_value;
	v2 = 0.0;
	TD::reset(initial_value);
}

// Linear Differentiator
void TD_LD::update(numeric_type ref) {
	/*
	 * Linear Differentiator.
	 * Parameter: lambda (bandwidth parameter, default 1)
	 */
	
	// Compute derivatives
	numeric_type dot_v1 = v2 - (v1 - ref) / lambda;
	numeric_type dot_v2 = -2.0 * (v1 - ref) / lambda;
	
	// Forward Euler integration
	v1 = v1 + dot_v1 * dt;
	v2 = v2 + dot_v2 * dt;
	
	// Outputs
	y = v1;
	yd = v2;
	ydd = 0.0;  // Second derivative not estimated
}

void TD_LD::reset(numeric_type initial_value) {
	v1 = initial_value;
	v2 = 0.0;
	TD::reset(initial_value);
}

TD_RED::TD_RED(numeric_type dt, numeric_type lambda1 = 1, numeric_type lambda2 = 1)
	: TD(dt)
{
	this->lambda1 = lambda1;
	this->lambda2 = lambda2;
}

// Robust Exact Differentiator
void TD_RED::update(numeric_type ref) {
	/*
	 * Robust Exact Differentiator.
	 * Parameters: lambda1, lambda2 (both default to 1)
	 */
	
	// Compute derivatives
	numeric_type err = v1 - ref;
	numeric_type dot_v1 = v2 - lambda1 * std::sqrt(std::abs(err)) * (err > 0 ? 1.0 : -1.0);
	numeric_type dot_v2 = -lambda2 * (err > 0 ? 1.0 : -1.0);
	
	// Forward Euler integration
	v1 = v1 + dot_v1 * dt;
	v2 = v2 + dot_v2 * dt;
	
	// Outputs
	y = v1;
	yd = v2;
	ydd = 0.0;  // Second derivative not estimated
}

void TD_RED::reset(numeric_type initial_value) {
	v1 = initial_value;
	v2 = 0.0;
	TD::reset(initial_value);
}

TD_INTD::TD_INTD(numeric_type dt, numeric_type alpha, numeric_type beta, numeric_type gamma, numeric_type r)
	:TD(dt)
{
	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;
	this->r = r;
}

// Improved Nonlinear Tracking Differentiator
void TD_INTD::update(numeric_type ref) {
	/*
	 * Improved Nonlinear Tracking Differentiator.
	 * Parameters: alpha, beta, gamma, r (defaults: 1, 1, 1, 1)
	 */
	
	// Parameter validation
	if (alpha < 0 || alpha > 1) {
		std::cerr << "Warning: alpha should be between 0 and 1 for INTD method, clamping" << std::endl;
		alpha = std::clamp<numeric_type>(alpha, 0.0, 1.0);
	}
	if (beta <= 0 || gamma <= 0 || r <= 0) {
		throw std::invalid_argument("beta, gamma, and r must be positive for INTD method");
	}
	
	// Compute derivatives
	numeric_type dot_v1 = v2;
	numeric_type dot_v2 = -r * r * std::tanh((beta * v1 - (1.0 - alpha) * ref) / gamma) 
					- r * v2;
	
	// Forward Euler integration
	v1 = v1 + dot_v1 * dt;
	v2 = v2 + dot_v2 * dt;
	
	// Outputs
	y = v1;
	yd = v2;
	ydd = 0.0;  // Second derivative not estimated
}

void TD_INTD::reset(numeric_type initial_value) {
	v1 = initial_value;
	v2 = 0.0;
	TD::reset(initial_value);
}