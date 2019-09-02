#pragma once

#include <stdint.h>
#include <vector>
#include <iostream>

class LongNumberIf
{
#if 1 // assertions
private:
	bool ValidZero() const;
	bool ValidDigits() const;
#endif // assertions
private:
	void TruncateZeros(std::vector<uint64_t>& digitsC) const;
	void AbsSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1, std::vector<uint64_t>& result) const;
	void AbsDifference(const std::vector<uint64_t>& minuend, const std::vector<uint64_t>& subtrahend, bool& resultSign, std::vector<uint64_t>& resultDigits) const;
private:
	LongNumberIf(const size_t& digitsN);
	LongNumberIf(const std::vector<uint64_t>& reversed);
	LongNumberIf DivShort(const uint32_t& div, uint32_t& remainder) const;

	std::vector<uint64_t> digits;
	bool negative = false;
public:
	LongNumberIf() = default;
	LongNumberIf(const std::vector<uint32_t>& externalDigits, bool isNegative);
	explicit LongNumberIf(const int64_t& number);
	explicit LongNumberIf(const uint32_t& number);
	LongNumberIf operator+(const LongNumberIf& other) const;
	LongNumberIf& operator+=(const LongNumberIf& other);
	bool operator==(const LongNumberIf& other) const;
	bool operator<(const LongNumberIf& other) const;
	friend std::ostream& operator<<(std::ostream& os, const LongNumberIf& number);
};