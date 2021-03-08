/* Author:  Mikhail Zhigun
   Email:   mikhail.zhigun@meteoswiss.ch  */

#include "Report.hpp"
#include <fstream>

namespace
{
	string trim(const string& str)
	{
		const size_t startIdx = str.find_first_not_of(" \t");
		if(startIdx == string::npos)
		{
			return {};
		}
		const size_t endIdx =  str.find_last_not_of(" \t");
		if(endIdx == string::npos)
		{
			return {};
		}
		return str.substr(startIdx, endIdx - startIdx + 1);
	}
}

extern "C"
{
	ecrad::Report* createReport(const char* filePath, int filePathLen)
	{
		const string path{filePath, static_cast<size_t>(filePathLen)};
		const string trimmedPath = trim(path);
		cout << __FUNCTION__ << " " << trimmedPath << endl;
		return new ecrad::Report{trimmedPath};
	}

	void deleteReport(ecrad::Report* report)
	{
		delete report;
	}

	void saveReport(ecrad::Report* report)
	{
		cout << __FUNCTION__ << " " << report->filePath() << endl;
		ofstream f{report->filePath()};
		if(f.is_open())
		{
			report->serialize(f);
		}
		else
		{
			cout << "Failed to open \"" <<  report->filePath() << "\"" << endl;
		}
		cout << __FUNCTION__ << " report saved" << endl;
	}

	void startReportStep(ecrad::Report* report, const char* name, int nameLen)
	{
		const string nameStr{name, static_cast<size_t>(nameLen)};
		cout << __FUNCTION__ << " " << nameStr << endl;
		report->startStep(trim(nameStr));
	}

	void finishCurrentReportStep(ecrad::Report* report, int result)
	{
		cout << __FUNCTION__  << " " << result << endl;
		cout << __FUNCTION__  << " " << report->steps().size() << endl;
		report->finishCurrentStep(result != 0);
	}

	void failCurrentReportStep(ecrad::Report* report, const char* name, int nameLen)
	{
		report->finishCurrentStep(false, trim(string{name, static_cast<size_t>(nameLen)}));
	}
}

namespace ecrad
{
	const Timestamp ReportStep::_programStartTS = chrono::high_resolution_clock::now();

	ReportStep::ReportStep(const string& name) :
		_name{ name },
		_startTS{ chrono::high_resolution_clock::now() }
	{}

	void ReportStep::finish(bool result, const std::string& errMsg)
	{
		_result = make_unique<bool>(result);
		_endtTS = make_unique<Timestamp>(chrono::high_resolution_clock::now());
		_errMsg = errMsg;
	}

	void ReportStep::serialize(ostream& ostrm) const
	{
		cout << __FUNCTION__ <<  endl;
		ostrm << "\t<step\n" <<
			"\t\tname=\"" << name() << "\"\n" <<
			"\t\tresult=\"" << (result()? "true" : "false") << "\"\n " <<
			"\t\tstart_timestamp_us=\"" << startTS().count() << "\"\n " <<
			"\t\tend_timestamp_us=\"" << endTS().count() << "\"";
		if(!result())
		{
			ostrm << "\n\t\terror_message=\"" << errorMessage().c_str() << "\"\n";
		}
		ostrm<<"/>\n";
	}

	inline bool ReportStep::finished() 
	{ 
		return static_cast<bool>(_result); 
	}

	inline chrono::microseconds ReportStep::startTS() const { return chrono::duration_cast<chrono::microseconds>(_startTS - _programStartTS); }

	inline chrono::microseconds ReportStep::endTS() const { return chrono::duration_cast<chrono::microseconds>((*_endtTS) - _programStartTS); }

	Report::Report(const string& filePath):
		_filePath{filePath}
	{}

	bool Report::startStep(const string& name)
	{
		if (!_steps.empty() && !_steps.back().finished())
		{
			return false;
		}
		_steps.emplace_back(name);
		cout << __FUNCTION__ << " " <<  steps().size() << endl;
	     return true;
	}

	bool Report::finishCurrentStep(bool result, const string& errorMessage)
	{
		if (_steps.empty() || _steps.back().finished())
		{
			return false;
		}
		_steps.back().finish(result, errorMessage);
		cout << __FUNCTION__ << " " <<  steps().size() << endl;
		return true;
	}

	void Report::serialize(ostream& ostrm) const
	{
		ostrm << "<?xml version=\"1.0\"?>\n" <<
			"<ecrad_driver_report>\n";
		cout << __FUNCTION__ << " " <<  steps().size() << endl;
		for (const auto& step : steps())
		{ step.serialize(ostrm); }
		ostrm << "</ecrad_driver_report>\n";
	}
}
