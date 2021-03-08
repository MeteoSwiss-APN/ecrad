/* Author:  Mikhail Zhigun
   Email:   mikhail.zhigun@meteoswiss.ch  */
#ifndef ECRAD_REPORT_HPP
#define ECRAD_REPORT_HPP

#include <iostream>
#include <string>
#include <chrono>
#include <memory>
#include <vector>

using namespace std;

namespace ecrad
{
	using Timestamp = chrono::time_point<chrono::high_resolution_clock>;

	class ReportStep
	{
	public:
		ReportStep(const string& name);
		void finish(bool result, const std::string& errMsg = std::string{});
		void serialize(ostream& ostrm) const;
		bool finished();
		const string& name() const { return _name; }
		const string& errorMessage() const { return _errMsg; }
		bool result() const { return *_result; }
		chrono::microseconds startTS() const;
		chrono::microseconds endTS() const;
	private:
		const string _name;
		std::string _errMsg;
		unique_ptr<bool> _result;
		Timestamp _startTS;
		unique_ptr<Timestamp> _endtTS;
		static const Timestamp _programStartTS;
	};

	class Report
	{
	public:
		Report(const string& filePath = string{});
		bool startStep(const string& name);
		bool finishCurrentStep(bool result, const string& errorMessage = std::string{});
		void serialize(ostream& ostrm) const;
		const vector<ReportStep>& steps() const { return _steps; }
		const string& filePath() const { return _filePath; }
	private:
		const string _filePath;
		vector<ReportStep> _steps;
	};
}

#endif // ECRAD_REPORT_HPP
