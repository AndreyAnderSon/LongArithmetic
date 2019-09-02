#pragma once

#include <stdint.h>
#include <vector>
#include <iostream>

class LongNumber
{
public:
	class DivideByZeroException : public std::exception
	{
	public:
		DivideByZeroException();
		~DivideByZeroException();
		virtual const char* what() const throw ();
	};
#if 1 // assertions
private:
	bool ValidZero() const;
	bool ValidDigits() const;
#endif // assertions
private:
	static void TruncateZeros(std::vector<uint64_t>& digitsC);
	void AbsSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1, std::vector<uint64_t>& result) const;
	static void AbsDifference(const std::vector<uint64_t>& minuend, const std::vector<uint64_t>& subtrahend, bool& resultSign, std::vector<uint64_t>& resultDigits);
	std::vector<uint64_t> GetSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1) const;
	static void LeftShift(std::vector<uint64_t>& numbers, uint8_t shift);
	static void RightShift(std::vector<uint64_t>& numbers, uint8_t shift);
private:
	LongNumber(const size_t& digits);
	LongNumber(const std::vector<uint64_t>& reversed);
	LongNumber(bool isNegative, const std::vector<uint64_t>& extDigits);
	LongNumber DivShort(const uint32_t& div, uint32_t& remainder) const;
		
	std::vector<uint64_t> digits;
	bool negative = false;
public:
	LongNumber() = default;
	LongNumber(const std::vector<uint32_t>& externalDigitsReversed, bool isNegative = false);
	explicit LongNumber(const int64_t& number);
	explicit LongNumber(const uint32_t& number);
	LongNumber(const std::string& number);

	LongNumber operator+(const LongNumber& other) const;
	LongNumber operator-(const LongNumber& other) const;
	LongNumber operator-() const;
	LongNumber Absolute() const;
	LongNumber operator*(const LongNumber& other) const;
	LongNumber QuotientRemainder(const LongNumber& denominator, LongNumber& remainder) const;
	LongNumber& operator+=(const LongNumber& other);
	LongNumber& operator=(const LongNumber& other);
	bool operator==(const LongNumber& other) const;
	bool operator!=(const LongNumber& other) const;
	bool operator<(const LongNumber& other) const;
	bool operator<=(const LongNumber& other) const;
	friend std::ostream& operator<<(std::ostream& os, const LongNumber& number);
#if 1 // gcd
public: 
	static LongNumber GCD(const LongNumber& ln1, const LongNumber& ln2);
private:
	static uint64_t CountZeroBitsRight(const std::vector<uint64_t>& number);
	static void LongRightShift(std::vector<uint64_t>& numbers, const uint64_t& shift);
	static void LongLeftShift(std::vector<uint64_t>& numbers, const uint64_t& shift);
#endif // gcd
};