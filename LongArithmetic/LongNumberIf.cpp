#include "LongNumberIf.h"
#include <algorithm>
#include <sstream>
#include <assert.h>

bool LongNumberIf::ValidZero() const
{
	return (digits.empty()) || (0 != digits.back());
}

bool LongNumberIf::ValidDigits() const
{
	bool result = true;
	for each (auto digit in digits) { result &= (0 == (digit & 0xffffffff00000000)); }
	return result;
}

void LongNumberIf::TruncateZeros(std::vector<uint64_t>& digitsC) const
{
	while (digitsC.size() > 0 && 0 == digitsC.back()) { digitsC.pop_back(); }
}

void LongNumberIf::AbsSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1, std::vector<uint64_t>& result) const
{
	assert(addend0.size() == addend1.size());
	assert(result.size() == addend0.size());
	assert(addend0.empty() || ((0 != addend0.back()) || (0 != addend1.back())));
	uint64_t carry = 0;
	for (size_t i = 0; i < addend0.size(); ++i)
	{
		uint64_t current = addend0[i] + addend1[i] + carry;
		carry = (current >> 32);
		assert(carry <= 1);
		result[i] = (current & 0x00000000ffffffff);
	}
	assert(carry <= 1);
	if (0 != carry) { result.push_back(carry); }
}

void LongNumberIf::AbsDifference(const std::vector<uint64_t>& minuend, const std::vector<uint64_t>& subtrahend, bool & resultSign, std::vector<uint64_t>& resultDigits) const
{
	assert(minuend.size() == subtrahend.size());
	assert(resultDigits.size() == subtrahend.size());
	assert(minuend.empty() || ((0 != minuend.back()) || (0 != subtrahend.back())));
	for (size_t i = 0; i < minuend.size(); ++i) { resultDigits[i] = minuend[i] - subtrahend[i]; }
	TruncateZeros(resultDigits);
	resultSign = !resultDigits.empty() && (0 != (resultDigits.back() & 0x8000000000000000));
	uint64_t add = static_cast<uint64_t>(resultSign);
	uint64_t xor (uint64_t(0) - add);
	uint64_t carry(0);
	for (auto digit = resultDigits.begin(); digit != resultDigits.end(); ++digit)
	{
		uint64_t current = (*digit ^ xor) + add + 0x0000000100000000 - carry;
		carry = uint64_t(1) - (current >> 32);
		assert(carry <= 1);
		*digit = (current & 0x00000000ffffffff);
	}
	assert(carry == 0);
	TruncateZeros(resultDigits);
}

LongNumberIf::LongNumberIf(const size_t & digitsN)
{
	digits.reserve(digitsN);
}

LongNumberIf::LongNumberIf(const std::vector<uint64_t>& reversed)
{
	digits.reserve(reversed.size());
	for (auto it = reversed.rbegin(); it != reversed.rend(); ++it) { digits.push_back(*it); }
	while (digits.size() > 0 && 0 == digits.back()) { digits.pop_back(); }
}

LongNumberIf LongNumberIf::DivShort(const uint32_t & div, uint32_t & remainder) const
{
	std::vector<uint64_t> reverse;
	reverse.reserve(digits.size());
	uint64_t rem = 0;
	for (auto it = digits.rbegin(); it != digits.rend(); ++it)
	{
		uint64_t current = (rem << 32) + *it;
		reverse.push_back(current / div);
		rem = current % div;
	}
	remainder = rem;
	return LongNumberIf(reverse);
}

LongNumberIf::LongNumberIf(const std::vector<uint32_t>& externalDigits, bool isNegative)
{
	digits.reserve(externalDigits.size());
	for (auto i = externalDigits.begin(); i != externalDigits.end(); ++i) { digits.push_back(*i); }
	TruncateZeros(digits);
}

LongNumberIf::LongNumberIf(const int64_t & number)
{
	negative = (number < 0);
	uint64_t abs(std::abs(number));
	digits.push_back(abs & 0x00000000ffffffff);
	digits.push_back(abs >> 32);
	TruncateZeros(digits);
}

LongNumberIf::LongNumberIf(const uint32_t & number)
{
	digits.clear();
	digits.push_back(number);
}

LongNumberIf LongNumberIf::operator+(const LongNumberIf & other) const
{
	LongNumberIf result;
	size_t longestNumberSize = std::max(digits.size(), other.digits.size());
	result.digits.reserve(longestNumberSize + 1);
	uint64_t carry = 0;
	for (size_t i = 0; i < longestNumberSize; ++i)
	{
		uint64_t current = digits[i] + other.digits[i] + carry;
		carry = (current >> 32);
		result.digits.push_back(current & 0x00000000ffffffff);
	}
	assert(false);
	//if (0 == carry){result.pus}
	return result;
}

LongNumberIf & LongNumberIf::operator+=(const LongNumberIf & other)
{
	assert(digits.empty() || (0 != digits.back()));
	assert(other.digits.empty() || (0 != other.digits.back()));
	size_t longestNumberSize = std::max(digits.size(), other.digits.size());
	std::vector<uint64_t> buffer = other.digits;
	digits.resize(longestNumberSize, 0);
	buffer.resize(longestNumberSize, 0);
	if (negative ^ other.negative)
	{
		bool order = (digits.empty() && !other.negative) ^ (!digits.empty() && negative);
		bool sign;
		AbsDifference(digits, buffer, sign, digits);
		negative = sign ^ order;
	}
	else
	{
		AbsSum(digits, buffer, digits);
	}
	return *this;
}

bool LongNumberIf::operator==(const LongNumberIf & other) const
{
	assert(ValidZero() && ValidDigits());
	assert(other.ValidZero() && other.ValidDigits());
	if (digits.empty() && other.digits.empty()) { return true; } // both are zeros
	else if ((digits.size() != other.digits.size()) || (negative ^ other.negative)) { return false; }
	else
	{
		bool result = true;
		for (size_t i = 0; i < digits.size(); ++i) { result &= (digits[i] == other.digits[i]); }
		return result;
	}
}

bool LongNumberIf::operator<(const LongNumberIf & other) const
{
	assert(ValidZero() && ValidDigits());
	assert(other.ValidZero() && other.ValidDigits());
	if (digits.empty() && other.digits.empty()) { return false; }
	else if (digits.size() != other.digits.size())
	{
		return
			(digits.size() < other.digits.size() && !other.negative) ||
			(digits.size() > other.digits.size() && negative);
	}
	else
	{
		bool lessAbs = false;
		bool equalAbs = true;
		for (int64_t i = digits.size() - 1; i >= 0; --i)
		{
			lessAbs |= (equalAbs && (digits[i] < other.digits[i]));
			equalAbs &= (digits[i] == other.digits[i]);
		}
		return (lessAbs && !other.negative) || (equalAbs && negative && !other.negative) || (!lessAbs && !equalAbs && negative);
	}
}

std::ostream & operator<<(std::ostream & os, const LongNumberIf & number)
{
	std::vector<uint32_t> tenDigits;
	LongNumberIf big = number;
	while (big.digits.size() != 0)
	{
		uint32_t cur;
		big = big.DivShort(1000000000, cur);
		tenDigits.push_back(cur);
	}
	auto it = tenDigits.rbegin();
	if (it == tenDigits.rend()) { os << 0; }
	else
	{
		os << *(it++);// << " ";
	}
	for (; it != tenDigits.rend(); ++it)
	{
		std::stringstream str;
		str << *it;
		std::string raw = str.str();
		std::string zeros(9 - raw.size(), '0');
		os << zeros << raw;// << " ";
	}
	return os;
}
