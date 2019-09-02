#include "LongNumber.h"
#include <algorithm>
#include <sstream>
#include <assert.h>
#include <stdexcept>

bool LongNumber::ValidZero() const
{
	return (digits.empty()) || (0 != digits.back());
}

bool LongNumber::ValidDigits() const
{
	bool result = true;
	//for each (auto digit in digits) { result &= (0 == (digit & 0xffffffff00000000)); }
	return result;
}

void LongNumber::TruncateZeros(std::vector<uint64_t>& digitsC)
{
	while (digitsC.size() > 0 && 0 == digitsC.back()) { digitsC.pop_back(); }
}

void LongNumber::AbsSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1, std::vector<uint64_t>& result) const
{
#if 1
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
#else
	assert(addend0.size() == addend1.size());
	assert(result.size() == addend0.size());
	assert(addend0.empty() || ((0 != addend0.back()) || (0 != addend1.back())));
	uint64_t carry = 0;
	for (size_t i = 0; i < addend0.size(); ++i)
	{
		uint64_t a0(addend0[i]);
		uint64_t a1(addend1[i]);
		uint64_t current = a0 + a1 + carry;
		result[i] = current;
		carry = uint64_t((current < a0) || (current < a1));
	}
	if (0 != carry) { result.push_back(carry); }
#endif
}

void LongNumber::AbsDifference(const std::vector<uint64_t>& minuend, const std::vector<uint64_t>& subtrahend, bool & resultSign, std::vector<uint64_t>& resultDigits)
{
	assert(minuend.size() == subtrahend.size());
	assert(resultDigits.size() == subtrahend.size());
	assert(minuend.empty() || ((0 != minuend.back()) || (0 != subtrahend.back())));
	for (size_t i = 0; i < minuend.size(); ++i) { resultDigits[i] = minuend[i] - subtrahend[i]; }
	TruncateZeros(resultDigits);
	resultSign = !resultDigits.empty() && (0 != (resultDigits.back() & 0x8000000000000000)); // subtrahend > minuend
	uint64_t add = static_cast<uint64_t>(resultSign);
	uint64_t xor(uint64_t(0) - add);
	uint64_t carry(0);
	for (auto digit = resultDigits.begin(); digit != resultDigits.end(); ++digit)
	{
		uint64_t current = (*digit ^ xor) + add + 0x0000'0001'0000'0000 - carry;
		carry = uint64_t(1) - (current >> 32);
		assert(carry <= 1);
		*digit = (current & 0x0000'0000'ffff'ffff);
	}
	assert(carry == 0);
	TruncateZeros(resultDigits);
}

std::vector<uint64_t> LongNumber::GetSum(const std::vector<uint64_t>& addend0, const std::vector<uint64_t>& addend1) const
{
	//std::vector<>
	return std::vector<uint64_t>();
}

void LongNumber::LeftShift(std::vector<uint64_t>& numbers, uint8_t shift)
{
	uint64_t carry = 0;
	for (auto current = numbers.begin(); current != numbers.end(); ++current)
	{
		*current = ((*current << shift) + carry);
		carry = (*current >> 32);
		*current &= 0x0000'0000'ffff'ffff;
	}
}

void LongNumber::RightShift(std::vector<uint64_t>& numbers, uint8_t shift)
{
	uint64_t carry = 0;
	for (auto current = numbers.rbegin(); current != numbers.rend(); ++current)
	{
		*current = ((*current << (32 - shift)) + carry);
		carry = (*current << 32);
		*current >>= 32;
	}
}

LongNumber::LongNumber(const size_t & digitsN)
{
	digits.reserve(digitsN);
}

LongNumber::LongNumber(const std::vector<uint64_t>& reversed)
{
	digits.reserve(reversed.size());
	for (auto it = reversed.rbegin(); it != reversed.rend(); ++it) { digits.push_back(*it); }
	TruncateZeros(digits);
}

LongNumber::LongNumber(bool isNegative, const std::vector<uint64_t>& extDigits) : negative(isNegative), digits(extDigits) {}

LongNumber LongNumber::DivShort(const uint32_t & div, uint32_t & remainder) const
{
	//LongNumber result(digits.size());
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
	return LongNumber(reverse);
}

LongNumber::LongNumber(const std::vector<uint32_t>& externalDigitsReversed, bool isNegative) : negative(isNegative)
{
	digits.reserve(externalDigitsReversed.size());
	for (auto i = externalDigitsReversed.rbegin(); i != externalDigitsReversed.rend(); ++i) { digits.push_back(*i); }
	TruncateZeros(digits);
}

LongNumber::LongNumber(const int64_t & number)
{
	negative = (number < 0);
	uint64_t abs(std::abs(number));
	digits.push_back(abs & 0x00000000ffffffff);
	digits.push_back(abs >> 32);
	TruncateZeros(digits);
}

LongNumber::LongNumber(const uint32_t & number)
{
	digits.clear();
	digits.push_back(number);
	// важно убрать нули, потому что в моей системе ноль задаётся только пустым вектором
	TruncateZeros(digits);
}

LongNumber::LongNumber(const std::string & number)
{
	LongNumber result;
	LongNumber base(uint32_t(1'000'000'000));
	bool neg = ('-' == number[0]);
	size_t start = neg ? 1 : 0;
	if (std::string::npos == number.find_first_not_of("0123456789", start))
	{
		size_t numberOfDecimalDigits = number.size() - start;
		size_t first = numberOfDecimalDigits % 9;
		result = result * base + LongNumber(uint32_t(std::stoi(number.substr(start, first))));
		start += first;
		for (size_t i = 0; i < numberOfDecimalDigits / 9; ++i)
		{
			result = result * base + LongNumber(uint32_t(std::stoi(number.substr(start, 9))));
			start += 9;
		}
	}
	*this = result;
	negative = neg;
}

LongNumber LongNumber::operator+(const LongNumber & other) const
{
	LongNumber result(negative, digits);
	result += other;
	return result;
}

LongNumber LongNumber::operator-(const LongNumber & other) const
{
	LongNumber result(!other.negative, other.digits);
	result += *this;
	return result;
}

LongNumber LongNumber::operator-() const
{
	return LongNumber(!negative, digits);
}

LongNumber LongNumber::Absolute() const
{
	return LongNumber(false, digits);
}

LongNumber LongNumber::operator*(const LongNumber & other) const
{
	assert(ValidZero() && ValidDigits());
	assert(other.ValidZero() && other.ValidDigits());
	std::vector<uint64_t> resultDigits(digits.size() + other.digits.size(), 0);
	const std::vector<uint64_t>* longer = (digits.size() > other.digits.size()) ? &digits : &other.digits;
	const std::vector<uint64_t>* shorter = (digits.size() > other.digits.size()) ? &other.digits : &digits;
	for (size_t j = 0; j < shorter->size(); ++j)
	{
		uint64_t carry = 0;
		for (size_t i = 0; i < longer->size(); ++i)
		{
			uint64_t current = longer->at(i) * shorter->at(j) + resultDigits[i + j] + carry;
			resultDigits[i + j] = (current & 0x0000'0000'ffff'ffff);
			carry = current >> 32;
		}
		assert(0 == resultDigits[j + longer->size()]);
		resultDigits[j + longer->size()] = carry;
	}
	TruncateZeros(resultDigits);
	return LongNumber(negative ^ other.negative, resultDigits);
}

LongNumber LongNumber::QuotientRemainder(const LongNumber & denominator, LongNumber & remainder) const
{
	assert(ValidZero() && ValidDigits());
	assert(denominator.ValidZero() && denominator.ValidDigits());
	// normalization
	if (denominator.digits.empty()) { throw DivideByZeroException(); }
	std::vector<uint64_t> denom = denominator.digits;
	std::vector<uint64_t> quotientReversed;
	std::vector<uint64_t> dividend = digits;
	dividend.resize(std::max(dividend.size(), denom.size()), 0);
	dividend.push_back(0);
	uint8_t shift = 32;
	for (uint64_t sample = denom.back(); sample > 0; sample >>= 1) { --shift; }
	assert(0x0000'0000'8000'0000 == ((denom.back() << shift) & 0xffff'ffff'8000'0000));
	LeftShift(dividend, shift);
	LeftShift(denom, shift);
	size_t N = dividend.size();
	size_t m = denom.size();
	quotientReversed.reserve(N - m);
	for (size_t i = 1; i <= (N - m); ++i)
	{
		uint64_t qEstimation = ((dividend[N - i] << 32) + dividend[N - i - 1]) / denom.back();
		assert(qEstimation <= uint64_t(0x0000'0000'ffff'ffff));
		dividend[N - i] += uint64_t(0xffff'ffff'0000'0000);
		uint64_t carry = 0;
		size_t begin = N - m - i;
		for (size_t j = 0; j < m; ++j)
		{
			dividend[begin + j] |= 0xffff'ffff'0000'0000;
			dividend[begin + j] -= (denom[j] * qEstimation + carry);
			carry = (dividend[begin + j] ^ 0xffff'ffff'ffff'ffff) >> 32;
			dividend[begin + j] &= 0x0000'0000'ffff'ffff;
		}
		dividend[N - i] -= carry;
		size_t limit = 2;
		while ((0xffff'ffff'0000'0000 != dividend[N - i]) && (limit > 0))
		{
			--limit;
			--qEstimation;
			uint64_t sCarry = 0;
			for (size_t j = 0; j < m; ++j)
			{
				dividend[begin + j] += (denom[j] + sCarry);
				sCarry = (dividend[begin + j] >> 32);
				dividend[begin + j] &= 0x0000'0000'ffff'ffff;
			}
			dividend[N - i] += sCarry;
		}
		assert(0xffff'ffff'0000'0000 == dividend[N - i]);
		dividend[N - i] = 0;
		quotientReversed.push_back(qEstimation);
	}
	LongNumber Q(quotientReversed);
	RightShift(dividend, shift);
	TruncateZeros(dividend);
	LongNumber R(false, dividend);
	//assert(R >= LongNumber(uint32_t(0)) && R < denominator);
	if (LongNumber(uint32_t(0)) != R && negative)
	{
		Q += LongNumber(uint32_t(1));
		R = LongNumber(false, denominator.digits) - R;
	}
	Q.negative = (negative ^ denominator.negative);
	remainder = R;
	return Q;
}

LongNumber & LongNumber::operator+=(const LongNumber & other)
{
#if 0 // realisation without if
	assert(ValidZero() && ValidDigits());
	assert(other.ValidZero() && other.ValidDigits());
	size_t numberOfDigits = std::max(digits.size(), other.digits.size());
	std::vector<uint64_t> a1 = digits;
	a1.resize(numberOfDigits, 0);
	std::vector<uint64_t> a2 = other.digits;
	a2.resize(numberOfDigits, 0);
	std::vector<uint64_t> resultDigits(numberOfDigits, 0);

	bool differentSigns = !digits.empty() && !other.digits.empty() && (negative ^ other.negative);
	{
		uint64_t add = static_cast<uint64_t>(differentSigns);
		uint64_t xor (uint64_t(0) - add);
		for (size_t i = 0; i < numberOfDigits; ++i) { resultDigits[i] = a1[i] + (a2[i] ^ xor) + add; }
		TruncateZeros(resultDigits);
	}
	bool negativeIntermediateResult = !resultDigits.empty() && (0 != (resultDigits.back() & 0x8000000000000000));
	bool resultSignNegative =
		(digits.empty() && other.negative) ||
		(other.digits.empty() && negative) ||
		(differentSigns && (negativeIntermediateResult ^ negative)) ||
		(!digits.empty() && !other.digits.empty() && other.negative && negative);
	uint64_t add = static_cast<uint64_t>(negativeIntermediateResult);
	uint64_t xor(uint64_t(0) - add);
	uint64_t carry(0);
	for (auto digit = resultDigits.begin(); digit != resultDigits.end(); ++digit)
	{
		uint64_t current = (*digit ^ xor) + add + 0x0000000100000000 - carry;
		carry = uint64_t(1) - (current >> 32);
		assert((1 == carry) || (0 == carry) || (0xffffffffffffffff == carry));
		*digit = (current & 0x00000000ffffffff);
	}
	if (0xffffffffffffffff == carry) { resultDigits.push_back(1); }
	TruncateZeros(resultDigits);
	negative = resultSignNegative;
	digits = resultDigits;
	return *this;
#else
	assert(ValidZero() && ValidDigits());
	assert(other.ValidZero() && other.ValidDigits());
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
#endif
}

LongNumber & LongNumber::operator=(const LongNumber & other)
{
	negative = other.negative;
	digits = other.digits;
	return *this;
}

bool LongNumber::operator==(const LongNumber & other) const
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

bool LongNumber::operator!=(const LongNumber & other) const
{
	return !(*this == other);
}

bool LongNumber::operator<(const LongNumber & other) const
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

bool LongNumber::operator<=(const LongNumber & other) const
{
	return (*this < other) || (*this == other);
}

LongNumber LongNumber::GCD(const LongNumber & ln1, const LongNumber & ln2)
{
	if (ln1.digits.empty()) { return ln2; }
	if (ln2.digits.empty()) { return ln1; }
	std::vector<uint64_t> buffer(1, 1);
	std::vector<uint64_t> M = ln1.digits;
	std::vector<uint64_t> N = ln2.digits;
	uint64_t mShift = CountZeroBitsRight(M);
	uint64_t nShift = CountZeroBitsRight(N);
	uint64_t commonShift = std::min(mShift, nShift);
	LongRightShift(M, mShift);
	LongRightShift(N, nShift);
	std::vector<uint64_t>* Mref = &M;
	std::vector<uint64_t>* Nref = &N;
	std::vector<uint64_t>* Bref = &buffer;
	std::vector<uint64_t>** biggest = &Bref;
	while (!(Mref->empty() || Nref->empty()))
	{
		LongRightShift(**biggest, CountZeroBitsRight(**biggest));
		size_t commonSize = std::max(Mref->size(), Nref->size());
		Mref->resize(commonSize, 0);
		Nref->resize(commonSize, 0);
		Bref->resize(commonSize, 0);
		bool N_bigger_M;
		AbsDifference(*Mref, *Nref, N_bigger_M, *Bref);
		biggest = N_bigger_M ? &Nref : &Mref;
		TruncateZeros(*reinterpret_cast<std::vector<uint64_t>*>(reinterpret_cast<intptr_t>(Mref) ^ reinterpret_cast<intptr_t>(Nref) ^ reinterpret_cast<intptr_t>(*biggest)));
		std::swap(*biggest, Bref);
	}
	assert(!(Mref->empty() && Nref->empty()));
	LongLeftShift(*Bref, commonShift);
	return LongNumber(false, *Bref);
}

uint64_t LongNumber::CountZeroBitsRight(const std::vector<uint64_t>& number)
{
	assert(!number.empty() && 0 != number.back());
	uint64_t zerosNumber = 0;
	auto digit = number.cbegin();
	for (; 0 == *digit; ++digit) { zerosNumber += 32; }
	uint64_t nzDigit = *digit;
	assert(nzDigit >> 32 == 0);
	while ((nzDigit & uint64_t(1)) == 0) 
	{
		nzDigit >>= 1;
		++zerosNumber;
	}
	return zerosNumber;
}

void LongNumber::LongRightShift(std::vector<uint64_t>& numbers, const uint64_t& shift)
{
	size_t start(shift >> 5);
	for (size_t i = start; i < numbers.size(); ++i) { numbers[i - start] = numbers[i]; }
	for (size_t i = numbers.size() - start; i < numbers.size(); ++i) { numbers[i] = 0; }
	RightShift(numbers, uint8_t(shift & 0x0000'0000'0000'001f));
	TruncateZeros(numbers);
}

void LongNumber::LongLeftShift(std::vector<uint64_t>& numbers, const uint64_t & shift)
{
	size_t add(shift >> 5);
	size_t initSize = numbers.size();
	numbers.resize(initSize + add + 1, 0);
	for (size_t i = 0; i < initSize; ++i) { numbers[i + add] = numbers[i]; }
	for (size_t i = 0; i < add; ++i) { numbers[i] = 0; }
	assert(0 == numbers.back());
	LeftShift(numbers, uint8_t(shift & 0x0000'0000'0000'001f));
	TruncateZeros(numbers);
}

std::ostream & operator<<(std::ostream & os, const LongNumber & number)
{
//	for (auto it = number.digits.rbegin(); it != number.digits.rend(); ++it)
//	{
//		os << *it << " ";
//	}
//	os << std::endl;
	std::vector<uint32_t> tenDigits;
	LongNumber big = number;
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

LongNumber::DivideByZeroException::DivideByZeroException()
{
}

LongNumber::DivideByZeroException::~DivideByZeroException()
{
}

const char * LongNumber::DivideByZeroException::what() const throw()
{
	return "LongNumber denominator has to be not 0";
}
