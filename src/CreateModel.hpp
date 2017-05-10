/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef CREATEMODEL_HPP_
#define CREATEMODEL_HPP_

#include "DrugDataReader.hpp"

#include "shannon_wang_puglisi_weber_bers_2004Cvode.hpp"
#include "fink_noble_giles_model_2008Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_endoCvode.hpp"
#include "ten_tusscher_model_2006_MCvode.hpp"
#include "noble_model_1998Cvode.hpp"
#include "ohara_rudy_2011_endoCvode.hpp"
#include "AbstractCvodeCell.hpp"
class CreateModel
{
public:
	unsigned mInterIndex;
	double mStartTime, mEndTime, mDt;
	boost::shared_ptr<AbstractCvodeCell> mpModel;
	boost::shared_ptr<RegularStimulus> mpStimulus;
	std::vector<double> mInterventionValues;
	std::vector<std::string> mDrugParameterNames, mIntervention;
	std::vector<double> mDrugConductanceFactors, mSteadyStateVariables, mOldSolution, mDruggedConductanceValues, mInterventionInitialState, mLimits;
	std::string mModelName, mInterventionType, mDrugName, mTraceFolder, mInterventionHardLimit;
	int mConcIndex;
    CreateModel(unsigned int inter_index, unsigned int model_index, boost::shared_ptr<RegularStimulus> p_stimulus, std::vector<double> intervention_values,
    		double start_time, double end_time, double dt, std::string drug_name, int conc_index, std::string trace_filename);
    CreateModel();
    void Drug (std::vector<double>conductance_factors, std::vector<std::string>channel_names);
    bool GetToSteadyState ();
    bool ResetAndCheck ();
    int Intervene (int c, unsigned conc_index, std::string drug_name);
    int InterveneAndCheck ();
    int DoTheIntervention (double intervention_value, std::string output_filename);
    std::string GetIdentifier();
    bool GetExpectedChange();
    bool IntervalBisection(bool result);
    void SetLimits(std::vector<double> coarse_threshold);
    double GetThreshold();
    std::vector<int> PossibleCombinations();
    std::string GetModelName();
    std::string GetInterventionName();


};

#endif
