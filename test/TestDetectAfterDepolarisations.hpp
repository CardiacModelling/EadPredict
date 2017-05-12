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

#ifndef TESTDETECTAFTERDEPOLARISATIONS_HPP_
#define TESTDETECTAFTERDEPOLARISATIONS_HPP_
#include <cxxtest/TestSuite.h>
#include "DetectAfterDepolarisations.hpp"
#include "ColumnDataReader.hpp"
#include "FileFinder.hpp"

class TestDetect : public CxxTest::TestSuite
{
public:
	void TestDetectAfterDepolarisations() throw(Exception)
	{
		/* Make sure this code is only run if CVODE is installed/enabled on the computer */
#ifdef CHASTE_CVODE


				double stim_start = 100;
				double stim_period = 3000;
				double lowlimit = 50;
				double highlimit = 100;

		FileFinder traces_folder("projects/BethM/test/TestTraces", RelativeTo::ChasteSourceRoot);
		TS_ASSERT_EQUALS(traces_folder.IsDir(), true);

		std::vector<FileFinder> traces = traces_folder.FindMatches("*.dat");

		for (unsigned i=0; i<traces.size(); i++)
		{
			std::cout << "Running for " << traces[i].GetLeafName() << std::endl;

			boost::shared_ptr<ColumnDataReader> p_reader(new ColumnDataReader(traces_folder,
					traces[i].GetLeafNameNoExtension()));
			std::vector<double> voltages = p_reader->GetValues("membrane_voltage");
			std::vector<double> simtimes = p_reader->GetValues("Time");
			double end_time = simtimes.back();
			std::cout << end_time;
			boost::shared_ptr<DetectAfterDepolarisations> p_AD(new DetectAfterDepolarisations);
			bool AD = p_AD->FindAD(voltages, simtimes, end_time, stim_start, stim_period, lowlimit, highlimit, "TestDetectAfterdepolarisations", true);
							if (traces[i].GetLeafName()[0] == 'N')
							{
								TS_ASSERT_EQUALS(AD,false);
							}
							else
							{
								TS_ASSERT_EQUALS(AD,true);
							}




		}

#else
		/* CVODE is not enable or installed*/
		std::cout << "Cvode is not enabled.\n";
#endif
	}
};

#endif /*TESTDETECTAFTERDEPOLARISATIONS_HPP_*/
