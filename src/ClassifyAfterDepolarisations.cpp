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


#include "ClassifyAfterDepolarisations.hpp"
#include <fstream>


ClassifyAfterDepolarisations::ClassifyAfterDepolarisations() {}

std::string ClassifyAfterDepolarisations::CausedAD(OdeSolution& solution, boost::shared_ptr<RegularStimulus> p_stimulus, std::string output_filename)
{
	/* get out voltage and time vectors */
	std::vector<double> v_Voltage = solution.GetAnyVariable("membrane_voltage");
	std::vector<double> v_Times = solution.rGetTimes();

	double end_time = v_Times.back();
	double stim_start = p_stimulus->GetStartTime();
	double stim_period = p_stimulus->GetPeriod();
	double TimeToIgnoreBeforeStimulus = 50;
	double TimeToIgnoreAfterStimulus = 100;

	return FindAD(v_Voltage, v_Times, end_time, stim_start, stim_period, TimeToIgnoreBeforeStimulus, TimeToIgnoreAfterStimulus, output_filename);
}

std::string ClassifyAfterDepolarisations::FindAD(std::vector<double> v_Voltage, std::vector<double> v_Times,
		double end_time, double stim_start, double stim_period, double TimeToIgnoreBeforeStimulus, double TimeToIgnoreAfterStimulus, std::string output_filename)
{

	std::vector<double> v_VoltageDifference;
	std::vector<double> v_TimesADIsOccurring;


	for (unsigned int a=0; a < v_Voltage.size()-1; a++)
	{
		/* list of voltage change at each time point */
		v_VoltageDifference.push_back(v_Voltage[a+1]-v_Voltage[a]);

		/* if dV/dt > 0.1mV/ms */
		if (v_VoltageDifference[a] > 0.001 )
		{
				/* store the time */
				v_TimesADIsOccurring.push_back(v_Times[a]);
		}
	}

	/* erase depolarisations provoked by stimulus */

	for (double stim_index = stim_start; stim_index < end_time; stim_index = stim_index + stim_period)
	{

		for (unsigned int ad_times_index = 0; ad_times_index < v_TimesADIsOccurring.size(); ad_times_index++)
		{

			if ((v_TimesADIsOccurring[ad_times_index] > stim_index -
					TimeToIgnoreBeforeStimulus && v_TimesADIsOccurring[ad_times_index] < stim_index + TimeToIgnoreAfterStimulus ) || v_TimesADIsOccurring[ad_times_index] < 100)
			{
				// delete
				v_TimesADIsOccurring.erase(v_TimesADIsOccurring.begin()+ad_times_index);
				// Go back one to make up for the deleted element
				ad_times_index--;
			}
		}

	}




	/* find times when the AD starts */
	std::vector<double> v_TimeDifference;
	std::vector<double> v_TimesADsBegin;

	/* iterate through list of voltage upswing times*/

	if (v_TimesADIsOccurring.size()>0)
	{

		for (unsigned int b=0; b < v_TimesADIsOccurring.size()-2; b++)
		{
			/* make a list of the differences */
			v_TimeDifference.push_back(v_TimesADIsOccurring[b+1]-v_TimesADIsOccurring[b]);
		}


		/* find minimum time difference */
		std::vector<double>::iterator step = std::min_element(v_TimeDifference.begin(),v_TimeDifference.end());
		double minsteps = *step;
		//std::cout << "MINSTEPS IS: "<< *step << "\n";

		for (unsigned int d=0; d< v_TimeDifference.size(); d++)
		{

			/* find all points in diff(adtimes) that change more than the minimum*/
			if (v_TimeDifference[d] > minsteps + 0.0001 )
			{
				/* store times */
				v_TimesADsBegin.push_back(v_TimesADIsOccurring[d+2]);
				for (unsigned int voltage_count = 0; voltage_count < v_Voltage.size() ; voltage_count++)
				{
					if (v_Times[voltage_count] == v_TimesADIsOccurring[d+2])
					{
					mVoltagesAtWhichADsStart.push_back(v_Voltage[voltage_count-1]);
					//std::cout << v_Voltage[voltage_count-1] << "\n";
					}
				}
			}
		}

	}


	std::vector<double> v_IgnoreTimes;
	int maximum_ADs = 0;
	double total_ADs = 0;
	std::vector<int> counting_ADs;
	for (double stim_index = stim_start; stim_index < end_time; stim_index = stim_index + stim_period)
	{

		int number_of_ADs = 0;

		for (unsigned int q=0; q<v_TimesADsBegin.size(); q++)
		{
			if (v_TimesADsBegin[q] > stim_index && v_TimesADsBegin[q] < stim_index+stim_period)
			{
				number_of_ADs += 1;
				total_ADs += 1;
			}
		}
		if (number_of_ADs > maximum_ADs)
		{
			maximum_ADs = number_of_ADs;
		}
		counting_ADs.push_back(number_of_ADs);
	}
	mADCount = counting_ADs;
	ADDensity = total_ADs/counting_ADs.size();
	mTotalAds = total_ADs;

	//std::cout << "Average number of ADs:" << ADDensity << "\n";

	std::ofstream outputfile;
	std::string file_location = getenv("CHASTE_TEST_OUTPUT");
	std::string file_name = file_location + "CADCheck/" + output_filename;
	outputfile.open(file_name.c_str());
	outputfile << "Time(s) membrane_voltage(mV) Ignore AD Repolarisation Recentness \n";

	bool recently_repolarised = false;
	bool after_AD = true;
	bool is_DAD = false;
	bool AD_occurred = false;
	bool repolarises = false;
	for (unsigned int d=0; d < v_Voltage.size(); d++)
	{
		outputfile << v_Times[d] << " " << v_Voltage[d] << " ";

		/*output 1 if AD, 0 if not */
		bool found = 0;
		for (unsigned int e=0; e < v_TimesADIsOccurring.size(); e++)
		{
			if (v_Times[d] == v_TimesADIsOccurring[e])
			{
				AD_occurred = true;
				after_AD = true;
				found = 1;
				if (recently_repolarised)
				{
					is_DAD = true;
					recently_repolarised = false;
				}
				break;
			}
		}
		outputfile << found << " ";
		bool ignore = 0;
		for (double stim_index = stim_start; stim_index < end_time; stim_index = stim_index + stim_period)
		{
			if (v_Times[d] == stim_index+TimeToIgnoreAfterStimulus)
			{
				recently_repolarised = false;
				after_AD = false;
			}

			if (v_Times[d] < stim_index+TimeToIgnoreAfterStimulus && v_Times[d] > stim_index-TimeToIgnoreBeforeStimulus)
			{
				ignore = 1;
				break;
			}
		}
		outputfile << ignore << " ";

		if (v_Voltage[d] < -60 && !after_AD)
		{
			recently_repolarised = true;
		}


		/* output 1 if repolarised, 0 if not */
		if (v_Voltage[d] < -60 && AD_occurred)
		{
			outputfile << "1 ";
			repolarises = true;
		}
		else
		{
			outputfile << "0 ";
		}

		outputfile << recently_repolarised << " ";
		outputfile << "\n";
	}
	outputfile.close();

	std::string AD_type = "E";

	if (is_DAD)
	{
		AD_type = "D";
	}

	if (!AD_occurred)
	{
		return "NNN";
	}
	else
	{
		if (repolarises)
		{
			if (maximum_ADs > 1)
				return AD_type + "MW";
			else
				return AD_type + "SW";
		}
		else
		{
			if (maximum_ADs > 1)
				return AD_type + "MO";
			else
				return AD_type + "SO";

		}
	}

	return "XXX";
}

int ClassifyAfterDepolarisations::GetCount()
{
	return mTotalAds;
}

std::vector<double> ClassifyAfterDepolarisations::GetStarts()
{
	return mVoltagesAtWhichADsStart;
}


#endif
