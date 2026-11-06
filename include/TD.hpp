#pragma once

#include <vector>
#include <string>
#include <map>

class TD {
public:
	typedef float numeric_type;
	// Output properties
	numeric_type y;      // Filtered output
	numeric_type yd;     // First derivative
	numeric_type ydd;    // Second derivative
	
	// Constructor
	TD(numeric_type dt);
	
	// Update method
	virtual void update(numeric_type ref) = 0;
	
	// Reset method
	virtual void reset(numeric_type initial_value = 0.0);
	
	// Set sample time
	void set_sample_time(numeric_type dt);

protected:
	numeric_type dt;
};

//Euler-based low passs filter with derivative estimation
class TD_Euler : public TD
{
public:
	//parameters
	numeric_type a; //filter coefficient
private:
	//state
	numeric_type prev_y, prev_yd;
public:
	TD_Euler(TD::numeric_type dt, TD::numeric_type a = 0.9);

	void update(numeric_type ref) final;

	void reset(numeric_type initial_value = 0.0) final;
};

//Time Optimal Differentiator
class TD_TOD : public TD
{
public:
	//parameters
	numeric_type r; //convergence rate
private:
	//state
	numeric_type v1, v2;
public:
	TD_TOD(numeric_type dt, numeric_type r = 1);

	void update(numeric_type ref) final;

	void reset(numeric_type initial_value = 0.0) final;
};

//Linear Differentiator
class TD_LD : public TD
{
public:
	//parameters
	numeric_type lambda; //bandwidth
private:
	//state
	numeric_type v1, v2;
public:
	TD_LD(numeric_type dt, numeric_type lambda = 1);

	void update(numeric_type ref) final;

	void reset(numeric_type initial_value = 0.0) final;
};

//Robust Exact Differentiator
class TD_RED : public TD
{
public:
	//parameters
	numeric_type lambda1, lambda2;
private:
	//state
	numeric_type v1, v2;
public:
	TD_RED(numeric_type dt, numeric_type lambda1 = 1, numeric_type lambda2 = 1);

	void update(numeric_type ref) final;

	void reset(numeric_type initial_value = 0.0) final;
};

//Improved nonlinear tracking differentiator
class TD_INTD : public TD
{
public:
	//parameters
	numeric_type alpha, beta, gamma, r;
private:
	//state
	numeric_type v1, v2;
public:
	TD_INTD(numeric_type dt, numeric_type alpha = 1.0, numeric_type beta = 1.0, numeric_type gamma = 1.0, numeric_type r = 1.0);

	void update(numeric_type ref) final;

	void reset(numeric_type initial_value = 0.0) final;
};