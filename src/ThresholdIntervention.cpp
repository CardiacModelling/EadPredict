/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifdef CHASTE_CVODE


#include "ThresholdIntervention.hpp"
#include <fstream>
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "DetectAfterDepolarisations.hpp"
#include "Debug.hpp"
#include "Debug.hpp"


ThresholdIntervention::ThresholdIntervention(){}

double ThresholdIntervention::FindAtUniformPace(boost::shared_ptr<AbstractCvodeCell> p_model,
		boost::shared_ptr<RegularStimulus> p_stimulus,  std::string protocol, int drug_index, int conc_index,
		double new_gNa_value, double new_gCaL_value, double new_gKr_value, double new_gpNa_value)
{
	// for running simulations
	double start_time = 0.0;
	double end_time = 10000.0;
	double dt = 0.1;

	bool result;
	double next_limit, starting_value, test_value;

	// initial values
	double threshold = 0;
	double limits[3] = {0,50,-2};
	int numsteps = 100;

	OdeSolution solution;
	std::string parameter;
	std::string modelname = p_model->GetSystemName();

	// reset to initial conditions and set currents to drugged state
	p_model->SetStateVariables(p_model->GetInitialConditions());
	p_model->SetParameter("membrane_fast_sodium_current_conductance",new_gNa_value);
	p_model->SetParameter("membrane_L_type_calcium_current_conductance",new_gCaL_value);
	p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance",new_gKr_value);
	try
	{
		p_model->SetParameter("membrane_persistent_sodium_current_conductance",new_gpNa_value);
	}
	catch(Exception &e)
	{
		std::cout << "No ipNa" << "\n";
	}

	// solve the model for a time period to check that resetting has been successful
	solution = p_model->Compute(start_time,end_time,dt);
	std::vector<double> old_solution = solution.GetAnyVariable("membrane_voltage");
	std::vector<double> new_solution;

	// set up intervention
	if (protocol=="iCaL")
	{
		parameter = "membrane_L_type_calcium_current_conductance";
		starting_value = p_model->GetParameter(parameter);
		limits[0] = 1;
	}
	if (protocol=="iNa")
	{
		parameter = "membrane_fast_sodium_current_shift_inactivation";
		starting_value = 1;
	}
	if (protocol == "iKr")
	{
		parameter = "membrane_rapid_delayed_rectifier_potassium_current_conductance";
		starting_value = p_model->GetParameter(parameter);
		limits[0] = -1;
		limits[1] = 1;
	}
	if (protocol == "Nai")
	{
		limits[0] = 1;
		limits[1] = 5;
	}

	for(int i=0; i<numsteps; i++)
	{
		if (protocol != "Nai")
		{
			//reset to initial state and set currents to drugged state
			p_model->SetStateVariables(p_model->GetInitialConditions());
			p_model->SetParameter("membrane_fast_sodium_current_conductance",new_gNa_value);
			p_model->SetParameter("membrane_L_type_calcium_current_conductance",new_gCaL_value);
			p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance",new_gKr_value);
			try
			{
				p_model->SetParameter("membrane_persistent_sodium_current_conductance",new_gpNa_value);
			}
			catch(Exception &e)
			{
				std::cout << "No ipNa" << "\n";
			}
			// reset intervention parameter
			p_model->SetParameter(parameter,starting_value);

			solution = p_model->Compute(start_time,end_time,dt);
			new_solution = solution.GetAnyVariable("membrane_voltage");

			// test that resetting has been successful
			if (i>0)
			{
				assert(old_solution == new_solution);
			}

			old_solution = new_solution;
			//reset to initial state again
			p_model->SetStateVariables(p_model->GetInitialConditions());
			p_model->SetParameter("membrane_fast_sodium_current_conductance",new_gNa_value);
			p_model->SetParameter("membrane_L_type_calcium_current_conductance",new_gCaL_value);
			p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance",new_gKr_value);
			try
			{
				p_model->SetParameter("membrane_persistent_sodium_current_conductance",new_gpNa_value);
			}
			catch(Exception &e)
			{
				std::cout << "No ipNa" << "\n";
			}
			/* change intervention parameter value */
			std::cout << "Iteration " << i << ": " << limits[0] << "\n";
			test_value = starting_value * limits[0];
			p_model->SetParameter(parameter,test_value);

			// for debugging
			PRINT_VARIABLE(p_model->GetParameter(parameter));
			PRINT_VARIABLE(p_model->GetParameter("membrane_fast_sodium_current_conductance"));
			PRINT_VARIABLE(p_model->GetParameter("membrane_L_type_calcium_current_conductance"));
			PRINT_VARIABLE(p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"));
			try
			{
				std::cout << "ipNa: " << p_model->GetParameter("membrane_persistent_sodium_current_conductance") << "\n";
			}
			catch(Exception &e)
			{
				std::cout << "No ipNa" << "\n";
			}
			/* run simulation */

			solution = p_model->Compute(start_time,end_time,dt);
		}

		else
		{
			// if modifying cytosolic sodium concentration
			//test_value = p_model[i]->GetStateVariable("cytosolic_sodium_concentration")*limits[0];
			p_model->SetStateVariable("cytosolic_sodium_concentration",limits[0]);
			boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
			p_regular_stim->SetPeriod(1000.0);
			OdeSolution early_solution = p_model->Compute(start_time,30000.0,dt);
			p_regular_stim->SetPeriod(3000.0);
			solution = p_model->Compute(30000.0,40000.0,dt);
		}

		// put together filename for output
		std::stringstream output_namestream;
		std::string output_filename;
		output_namestream << "CvodeCells" << protocol << modelname << "_s_" <<	drug_index << "_" << conc_index << "_" << limits[0];
		output_filename = output_namestream.str();
		/* is there an AD? */
		DetectAfterDepolarisations AD_test;
		result = AD_test.CausedAD(solution, p_stimulus, output_filename, false);

		// based on whether there is an AD, pick a new value to test

		if (limits[2]==-2) //first iteration
		{
			// limits = {lower_limit, upper_limit, -2}
			double hold_upper_limit = limits[1];
			limits[1] = limits[0];
			limits[0] = hold_upper_limit;
			limits[2] = -3;
			// limits = {upper_limit, lower_limit, -3}
			if (result)
			{
				limits[1]+=1;
			}
		}
		else if (limits[2]==-3) // second iteration
		{
			// limits = {upper_limit, lower_limit, -3}
			limits[2] = limits[0];
			limits[0] = limits[1]+((limits[2]-limits[1])/2);
			// limits = {middle_number, lower_limit, upper_limit}
		}
		else
		{
			if (result) // if there is an AD
			{
				//calculate middle number, get rid of old upper limit
				// limits = {middle_number, lower_limit, upper_limit}
				next_limit = limits[1]+((limits[0]-limits[1])/2);
				limits[2] = limits[0];
				limits[0] = next_limit;
				threshold = limits[0];
			}
			else
			{
				// calculate middle number, get rid of old lower limit
				next_limit = limits[0]+((limits[2]-limits[0])/2);
				limits[1] = limits[0];
				limits[0] = next_limit;

			}
			if ((limits[2]-limits[1])< 0.0001) // if thresholds are converging
			{
				break;
			}

		}
	}


	return threshold;
}


#endif
