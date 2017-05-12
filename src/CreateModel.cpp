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

#include "DrugDataReader.hpp"
#include "SteadyStateRunner.hpp"
#include "CreateModel.hpp"
#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "fink_noble_giles_model_2008Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_endoCvode.hpp"
#include "ten_tusscher_model_2006_MCvode.hpp"
#include "noble_model_1998Cvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "grandi_pasqualini_bers_2010_ssCvode.hpp"
#include "AbstractCvodeCell.hpp"
#include "DetectAfterDepolarisations.hpp"
#include "ClassifyAfterDepolarisations.hpp"

#include "Timer.hpp"

CreateModel::CreateModel(unsigned int inter_index, unsigned int model_index,
		boost::shared_ptr<RegularStimulus> p_stimulus,
		std::vector<double> intervention_values,
		double start_time, double end_time, double dt, std::string drug_name, int conc_index,
		std::string trace_folder)
:  mInterIndex(inter_index),
   mStartTime(start_time),
   mEndTime(end_time),
   mDt(dt),
   mpStimulus(p_stimulus),
   mInterventionValues(intervention_values),
   mDrugName(drug_name),
   mTraceFolder(trace_folder),
   mConcIndex(conc_index)
{
	boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

	// pick model
	boost::shared_ptr<AbstractCvodeCell> p_model;
	if (model_index == 1)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellshannon_wang_puglisi_weber_bers_2004FromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}
	else if (model_index == 2)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellfink_noble_giles_model_2008FromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}
	else if (model_index == 3)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}
	else if (model_index == 4)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_endoFromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}
	else if (model_index == 5)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellten_tusscher_model_2006_MFromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}

	else if (model_index == 6)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellohara_rudy_2011_endoFromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;
	}
	else if (model_index == 7)
	{
		boost::shared_ptr<AbstractCvodeCell> p_new_model(new Cellgrandi_pasqualini_bers_2010_ssFromCellMLCvode(p_solver, mpStimulus));
		p_model = p_new_model;

	}

	mpModel = p_model;
	mModelName = mpModel->GetSystemName();

	// set up intervention
	if (mInterIndex == 1)
	{
		mInterventionType = "StateVariable";
		mIntervention.push_back("cytosolic_sodium_concentration");
		mInterventionHardLimit = "AboveZero";

	}
	else if (mInterIndex ==  2)
	{
		mInterventionType = "Parameter";
		mIntervention.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
		mInterventionHardLimit = "BelowOne";

	}
	else if (mInterIndex == 3)
	{
		mInterventionType = "Parameter";
		mIntervention.push_back("membrane_fast_sodium_current_shift_inactivation");
		mInterventionHardLimit = "AboveZero";

	}
	else if (mInterIndex == 4)
	{
		mInterventionType = "Parameter";
		mIntervention.push_back("membrane_L_type_calcium_current_conductance");
		mInterventionHardLimit = "AboveZero";

	}
	else if (mInterIndex == 5)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_L_type_calcium_current_conductance");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "AboveZero";
	}

	else if (mInterIndex == 6)
	{
		mInterventionType = "Parameter";
		mIntervention.push_back("membrane_persistent_sodium_current_conductance");
	}

	else if (mInterIndex == 7)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "BelowOne";
	}
	else if (mInterIndex == 8)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_fast_sodium_current_shift_inactivation");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "AboveZero";
	}
	else if (mInterIndex == 9)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_L_type_calcium_current_conductance");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "AboveZero";
	}
	else if (mInterIndex == 10)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "BelowOne";
	}
	else if (mInterIndex == 11)
	{
		mInterventionType = "TwoParameters";
		mIntervention.push_back("membrane_fast_sodium_current_shift_inactivation");
		mIntervention.push_back("extracellular_potassium_concentration");
		mInterventionHardLimit = "AboveZero";
	}
	else if (mInterIndex ==  12)
	{
		mInterventionType = "Parameter";
		mIntervention.push_back("membrane_fast_sodium_current_reduced_inactivation");
		mInterventionHardLimit = "AboveZero";

	}

}

CreateModel::CreateModel()
{
	std::cout << "Empty model created\n";
}

#endif

void CreateModel::Drug (std::vector<double>conductance_factors, std::vector<std::string> channel_names)
{
	// take conductance factors and apply them to the model
	mDrugConductanceFactors = conductance_factors;
	mDrugParameterNames = channel_names;
	std::vector<double> original_values;
	for (unsigned int i= 0; i<mDrugParameterNames.size(); i++)
	{
		try
		{
			original_values.push_back(mpModel->GetParameter(mDrugParameterNames[i]));
		}
		catch(Exception &e)
		{
			std::cout << "No " << mDrugParameterNames[i] << "\n";
		}
	}
	for (unsigned int j= 0; j<mDrugParameterNames.size(); j++)
	{
		try
		{
			mpModel->SetParameter(mDrugParameterNames[j],original_values[j]*mDrugConductanceFactors[j]);
			mDruggedConductanceValues.push_back(original_values[j]*mDrugConductanceFactors[j]);
		}
		catch(Exception &e)
		{
			std::cout << "\n";
		}
	}

	// get initial intervention values
	if (mInterIndex == 1)
	{
		mInterventionInitialState.push_back(1);
		std::cout << "Original value is: " << mpModel->GetStateVariable(mIntervention[0]) << "\n";
	}
	else if (mInterIndex ==  2)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
	}
	else if (mInterIndex == 3)
	{
		this->mInterventionInitialState.push_back(0);
	}
	else if (mInterIndex == 4)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));

	}
	else if (mInterIndex == 5)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));
	}
	else if (mInterIndex == 6)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
	}
	else if (mInterIndex == 7)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));
	}
	else if (mInterIndex == 8)
	{
		mInterventionInitialState.push_back(0);
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));

	}
	else if (mInterIndex == 9)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));
	}
	else if (mInterIndex == 10)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));
	}
	else if (mInterIndex == 11)
	{
		mInterventionInitialState.push_back(0);
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[1]));
	}

	else if (mInterIndex == 12)
	{
		mInterventionInitialState.push_back(mpModel->GetParameter(mIntervention[0]));
	}

}

bool CreateModel::GetToSteadyState()
{
	SteadyStateRunner steady_runner(mpModel);
	bool result = steady_runner.RunToSteadyState();
	
	// Reset intervention variable
	if (mInterventionType == "StateVariable")
	{
		mpModel->SetStateVariable(mIntervention[0],mInterventionInitialState[0]);
	}
	else if (mInterventionType == "Parameter")
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]);
	}
	else if (mInterventionType == "TwoParameters")
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]);

	}
	// Get all the variables out
	mSteadyStateVariables = mpModel->GetStdVecStateVariables();
	// For checking resetting later on
	OdeSolution solution = mpModel->Compute(mStartTime,mEndTime,mDt);
	mOldSolution = solution.GetAnyVariable("membrane_voltage");
	return result;
}

bool CreateModel::ResetAndCheck()
{

	// apply drug effects to model
	for (unsigned int i = 0 ; i< mDruggedConductanceValues.size(); i++)
	{
		mpModel->SetParameter(mDrugParameterNames[i], mDruggedConductanceValues[i]);

	}
	// reset state variables, including intervention state variable
	mpModel->SetStateVariables(mSteadyStateVariables);

	// reset intervention parameters if necessary
	if (mInterventionType == "Parameter")
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]);
	}
	else if (mInterventionType == "TwoParameters")
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]);
	}
	bool result = true;
	return result;
}

int CreateModel::Intervene(int c, unsigned conc_index, std::string drug_name)
{

	std::string output_filename;
	std::stringstream output_namestream;
	output_namestream << mTraceFolder << "/Sequential/" << mDrugName << mConcIndex << GetIdentifier() << c;
	output_filename = output_namestream.str();
	double intervention_value = mInterventionValues[c];
	return DoTheIntervention(intervention_value, output_filename);

}

int CreateModel::DoTheIntervention(double intervention_value, std::string output_filename)
{
	// apply the intervention to the model
	if (mInterIndex == 1)
	{
		mpModel->SetStateVariable(mIntervention[0],intervention_value);
		// check for a DELAYED afterdepolarisation
		mpStimulus->SetPeriod(1000.0);
		mpModel->Solve(0,30000,mpStimulus->GetDuration()); // Stimulus duration is the maximum time step to use
		mpStimulus->SetPeriod(3000.0);
		OdeSolution solution = mpModel->Solve(30000+mStartTime,30000+mEndTime,mpStimulus->GetDuration(),mDt);
		ClassifyAfterDepolarisations AD_classify;
		std::string result = AD_classify.CausedAD(solution, mpStimulus, output_filename);
		if (result[0] == 'D')
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else if (mInterIndex ==  2)
	{
		mpModel->SetParameter(mIntervention[0],(mInterventionInitialState[0])*intervention_value);
	}
	else if (mInterIndex == 3)
	{
		mpModel->SetParameter(mIntervention[0],intervention_value);
	}
	else if (mInterIndex == 4)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
	}

	else if (mInterIndex == 5)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*2);
	}
	else if (mInterIndex == 6)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
	}
	else if (mInterIndex == 7)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*2);
	}
	else if (mInterIndex == 8)
	{
		mpModel->SetParameter(mIntervention[0],intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*2);
	}
	else if (mInterIndex == 9)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*0.75);
	}
	else if (mInterIndex == 10)
	{
		mpModel->SetParameter(mIntervention[0],mInterventionInitialState[0]*intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*0.75);
	}
	else if (mInterIndex == 11)
	{
		mpModel->SetParameter(mIntervention[0],intervention_value);
		mpModel->SetParameter(mIntervention[1],mInterventionInitialState[1]*0.75);
	}
	else if (mInterIndex == 12)
	{
		mpModel->SetParameter(mIntervention[0],intervention_value);
	}
	// check for an afterdepolarisation
	OdeSolution solution ;

	try
	{
	    Timer::Reset();
		solution = mpModel->Solve(mStartTime,mEndTime,mpStimulus->GetDuration(), mDt);
		std::cout << mEndTime - mStartTime << " ms solve took " << Timer::GetElapsedTime() << "s.\n"<< std::flush;
	}
	catch (Exception &e)
	{
		std::cout << e.GetShortMessage() << std::endl;
		return -2;
	}

	Timer::Reset();
	DetectAfterDepolarisations AD_test;
	int result = AD_test.CausedAD(solution, mpStimulus, output_filename, true);
    std::cout << "AD detection took " << Timer::GetElapsedTime() << "s.\n"<< std::flush;

	return result;
}

std::string CreateModel::GetIdentifier()
{
	// return the name of the model and the intervention
	if (mInterIndex == 5)
	{
		return (mModelName + mIntervention[0] + "_double_" + mIntervention[1]);
	}
	else if (mInterIndex == 7)
	{
		return (mModelName + mIntervention[0] + "_double_" + mIntervention[1]);
	}
	else if (mInterIndex == 8)
	{
		return (mModelName + mIntervention[0] + "_double_" + mIntervention[1]);
	}
	else if (mInterIndex == 9)
	{
		return (mModelName + mIntervention[0] + "_three_quarters_" + mIntervention[1]);
	}
	else if (mInterIndex == 10)
	{
		return (mModelName + mIntervention[0] + "_three_quarters_" + mIntervention[1]);
	}

	else if (mInterIndex == 11)
	{
		return (mModelName + mIntervention[0] + "_three_quarters_" + mIntervention[1]);
	}

	return (mModelName + mIntervention[0]);
}


bool CreateModel::GetExpectedChange()
{
	// Returns true if the intervention is expected to go from no AD to an AD as the intervention increases (i.e. IKr conductance)
	// Returns false otherwise
	if (mInterIndex == 2 || mInterIndex == 7 || mInterIndex == 10 )
	{
		return false;
	}
	return true;
}

void CreateModel::SetLimits(std::vector<double> coarse_threshold)

{
	assert(coarse_threshold.size()==2);
	mLimits = coarse_threshold;
	mLimits.push_back(-2);
}

int CreateModel::InterveneAndCheck()
{
	std::stringstream output_namestream;
	std::string output_filename;
	output_namestream << mTraceFolder << "/Bisection/"  << mDrugName << mConcIndex << GetIdentifier() << mLimits[0];
	output_filename = output_namestream.str();
	double intervention_value = mLimits[0];
	return DoTheIntervention(intervention_value, output_filename);
}

bool CreateModel::IntervalBisection(bool result)
{
	// put together filename for AD checking
	std::stringstream output_namestream;
	std::string output_filename;
	output_namestream << mTraceFolder << "/Bisection/" << mDrugName << mConcIndex << GetIdentifier();	/* is there an AD? */
	output_filename = output_namestream.str();
	std::cout << "Intervention value:  " << mLimits[0] << "\n2:  " << mLimits[1] << "\n3:  " << mLimits[2] << "\n";

	bool expected_result = GetExpectedChange();

	// initial values
	if (mLimits[2]==-2) //first iteration
	{
		// mLimits = {lower_limit, upper_limit, -2}
		if (result == expected_result)
		{
			std::cout << "Unexpected result at lower limit\n";
			mLimits[0] = mLimits[0] - 1;
			// mLimits = {new_lower_limit, upper_limit, -2}

		}
		else
		{
			double hold_upper_limit = mLimits[1];
			mLimits[1] = mLimits[0];
			mLimits[0] = hold_upper_limit;
			mLimits[2] = -3;
			// mLimits = {upper_limit, lower_limit, -3}

		}
	}
	else if (mLimits[2]==-3) // second iteration
	{
		// mLimits = {upper_limit, lower_limit, -3}
		if (result != expected_result)
		{
			std::cout << "Unexpected result at higher limit\n";
			mLimits[0] = mLimits[0] + 1;
			// mLimits = {new_upper_limit, lower_limit, -3}

		}
		else
		{
			// mLimits = {upper_limit, lower_limit, -3}
			mLimits[2] = mLimits[0];
			mLimits[0] = mLimits[1]+((mLimits[2]-mLimits[1])/2);
			// mLimits = {middle_number, lower_limit, upper_limit}

		}
	}
	else // all other iterations
	{
		if (result == expected_result) // if there is an AD
		{
			//calculate middle number, get rid of old upper limit
			// mLimits = {middle_number, lower_limit, upper_limit}
			double next_limit = mLimits[1]+((mLimits[0]-mLimits[1])/2);
			mLimits[2] = mLimits[0];
			mLimits[0] = next_limit;
		}
		else
		{
			// calculate middle number, get rid of old lower limit
			double next_limit = mLimits[0]+((mLimits[2]-mLimits[0])/2);
			mLimits[1] = mLimits[0];
			mLimits[0] = next_limit;

		}
		if (std::abs(mLimits[2]-mLimits[1])< 0.0001) // if thresholds are converging
		{
			return true;
		}

	}
	return false;
}

double CreateModel::GetThreshold()
{
	return mLimits[1];
}

std::vector<int> CreateModel::PossibleCombinations()
{
	// the idea is to return the total number of models and total number of interventions so that
	// TestCompletenessOfData.hpp can work out how much data we don't have yet
	// currently it's hard-coded
	std::vector<int> inputs;
	inputs.push_back(5);
	inputs.push_back(12);
	//inputs[0] = 5;
	//inputs[1] = 12;
	return inputs;
}

std::string CreateModel::GetModelName()
{
	return mModelName;
}

std::string CreateModel::GetInterventionName()
{
	// the idea is to return the name of the intervention we are using
	// currently this is not implemented
	return "hello";
}

