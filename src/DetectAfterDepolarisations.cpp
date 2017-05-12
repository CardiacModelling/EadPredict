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

#ifdef CHASTE_CVODE


#include "DetectAfterDepolarisations.hpp"
#include <fstream>


DetectAfterDepolarisations::DetectAfterDepolarisations() {}

bool DetectAfterDepolarisations::CausedAD(OdeSolution& solution, boost::shared_ptr<RegularStimulus> p_stimulus, std::string output_filename, bool print_trace)
{
	/* get out voltage and time vectors */
	std::vector<double> v_Voltage = solution.GetAnyVariable("membrane_voltage");
	std::vector<double> v_Times = solution.rGetTimes();

	double end_time = v_Times.back();
	double stim_start = p_stimulus->GetStartTime();
	double stim_period = p_stimulus->GetPeriod();
	double TimeToIgnoreBeforeStimulus = 50;
	double TimeToIgnoreAfterStimulus = 100;

	return FindAD(v_Voltage, v_Times, end_time, stim_start, stim_period, TimeToIgnoreBeforeStimulus, TimeToIgnoreAfterStimulus, output_filename, print_trace);
}

bool DetectAfterDepolarisations::FindAD(std::vector<double> v_Voltage, std::vector<double> v_Times,
		double end_time, double stim_start, double stim_period, double TimeToIgnoreBeforeStimulus, double TimeToIgnoreAfterStimulus, std::string output_filename, bool print_trace)
{

	try
	{

		std::vector<double> v_VoltageDifference;
		std::vector<double> v_TimesADIsOccurring;


		for (unsigned int a=0; a < v_Voltage.size()-1; a++)
		{
			/* list of voltage change at each time point */
			v_VoltageDifference.push_back(v_Voltage[a+1]-v_Voltage[a]);

			/* if dV/dt > 1mV/ms */
			if (v_VoltageDifference[a] > 0.001)
				/* store the time */
				v_TimesADIsOccurring.push_back(v_Times[a]);
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
					/*                std::cout << adtimes[d+2] << "\n";*/
				}
			}

		}

		if (print_trace)
		{
			std::ofstream outputfile;

			std::string file_location = getenv("CHASTE_TEST_OUTPUT");
			std::string file_name = file_location + output_filename;
			outputfile.open(file_name.c_str());
			outputfile << "Time(s) membrane_voltage(mV) AD\n";
			for (unsigned int d=0; d < v_Voltage.size(); d++)
			{
				outputfile << v_Times[d] << " " << v_Voltage[d] << " ";

				/*output 1 if AD, 0 if not */
				bool found = 0;
				for (unsigned int e=0; e < v_TimesADIsOccurring.size(); e++)
				{
					if (v_Times[d] == v_TimesADIsOccurring[e])
					{
						found = 1;
						break;
					}
				}
				outputfile << found;

				outputfile << "\n";
			}
			outputfile.close();
		}
		if (v_TimesADsBegin.size()>0)
		{
			return 1;
		}

	}
	catch(Exception& e)
	{
		std::cout << "DetectAfterdepolarisations Failed";
		return 0;
	}


	return 0;
}
#endif
