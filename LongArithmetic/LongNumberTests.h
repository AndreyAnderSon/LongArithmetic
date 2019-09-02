#pragma once
#include <vector>
#include <random>
class LongNumber;

class LongNumberTests
{
private:
	class LongNumberComparatorTest
	{
	private:
		std::vector<std::vector<uint32_t>> GetVariants(const std::vector<std::vector<uint32_t>>& base) const;
		int GetNumber(const std::vector<uint32_t>& digits, bool negative) const;
		bool TestForThisCouple(bool sign0, const std::vector<uint32_t>& digits0, bool sign1, const std::vector<uint32_t>& digits1) const;
	public:
		bool Run() const;
	};
	class LongNumberAdditionTest
	{
	private:
		bool CheckPair(const int64_t& first, const int64_t& second) const;
		std::vector<int64_t> GetRandomNumbers(const size_t& count) const;
	public:
		bool Run() const;
	};
	class LongNumberMultiplicationTest
	{
	public:
		bool Run() const;
	};
	class LongNumberDivisionTest
	{
	private:
		std::mt19937_64 absState;
		bool CheckPair(const LongNumber& first, const LongNumber& second) const;
		bool CheckPairSigns(const LongNumber& first, const LongNumber& second) const;
		LongNumber GetRandomNumber();
	public:
		LongNumberDivisionTest();
		bool Run();
	};
	class LongNumberGCDTest
	{
	private:
		std::mt19937_64 absState;
		LongNumber GCD(const LongNumber& first, const LongNumber& second) const;
		bool CheckTrio(const LongNumber& first, const LongNumber& second, const LongNumber& common) const;
		LongNumber GetRandomNumber();
	public:
		LongNumberGCDTest();
		bool Run();
	};
	class PerformanceCompare
	{
	private:
		std::vector<int64_t> GetRandomNumbers(const size_t& count) const;
	public:
		void Run() const;
	};
public:
	static bool TestLongNumbers();
};