#include "LongNumberTests.h"
#include "LongNumber.h"
#include "LongNumberIf.h"
#include <random>
#include <limits>

std::vector<std::vector<uint32_t>> LongNumberTests::LongNumberComparatorTest::GetVariants(const std::vector<std::vector<uint32_t>>& base) const
{
	std::vector<std::vector<uint32_t>> result;
	result.reserve(base.size() * 4);
	for each (auto seed in base)
	{
		for (int i = 0; i < 4; ++i)
		{
			result.push_back(seed);
			result.back().push_back(i);
		}
	}
	return result;
}

int LongNumberTests::LongNumberComparatorTest::GetNumber(const std::vector<uint32_t>& digits, bool negative) const
{
	int result = 0;
	int pos = 1;
	for (auto i = digits.rbegin(); i != digits.rend(); ++i)
	{
		result += (pos * *i);
		pos *= 4;
	}
	return negative ? -result : result;
}

bool LongNumberTests::LongNumberComparatorTest::TestForThisCouple(bool sign0, const std::vector<uint32_t>& digits0, bool sign1, const std::vector<uint32_t>& digits1) const
{
	int n0 = GetNumber(digits0, sign0);
	int n1 = GetNumber(digits1, sign1);
	LongNumber ln0(digits0, sign0);
	LongNumber ln1(digits1, sign1);
	return ((n0 < n1) == (ln0 < ln1)) && ((n0 == n1) == (ln0 == ln1));
}

bool LongNumberTests::LongNumberComparatorTest::Run() const
{
	bool result = true;
	std::vector<std::vector<uint32_t>> numbers;
	numbers.push_back(std::vector<uint32_t>());
	for (size_t digitNumber = 0; digitNumber < 3; ++digitNumber) { numbers = GetVariants(numbers); }
	for each (auto first in numbers)
	{
		for each (auto second in numbers)
		{
			result &= TestForThisCouple(false, first, false, second);
			result &= TestForThisCouple(false, first, true, second);
			result &= TestForThisCouple(true, first, false, second);
			result &= TestForThisCouple(true, first, true, second);
		}
	}
	return result;
}

bool LongNumberTests::TestLongNumbers()
{
	bool result = true;
	{
		LongNumberComparatorTest comparatorTest;
		result &= comparatorTest.Run();
		std::cout << "LongNumberComparatorTest: " << result << std::endl;
	}
	{
		LongNumberAdditionTest additionTest;
		result &= additionTest.Run();
		std::cout << "LongNumberAdditionTest: " << result << std::endl;
	}
	{
		LongNumberMultiplicationTest multiplicationTest;
		result &= multiplicationTest.Run();
		std::cout << "LongNumberMultiplicationTest: " << result << std::endl;
	}
	{
		LongNumberDivisionTest divisionTest;
		result &= divisionTest.Run();
		std::cout << "LongNumberDivisionTest: " << result << std::endl;
	}
	{
		LongNumberGCDTest gcdTest;
		result &= gcdTest.Run();
		std::cout << "LongNumberGCDTest: " << result << std::endl;
	}
	//PerformanceCompare perf;
	//perf.Run();

	return result;
}

bool LongNumberTests::LongNumberAdditionTest::CheckPair(const int64_t & first, const int64_t & second) const
{
	LongNumber ln(first);
	ln += LongNumber(second);
	return ln == LongNumber(first + second);
}

std::vector<int64_t> LongNumberTests::LongNumberAdditionTest::GetRandomNumbers(const size_t& count) const
{
	std::random_device absSeed;
	std::mt19937_64 absState(absSeed());
	std::uniform_int_distribution<uint64_t> absDistr(uint64_t(0), uint64_t(0x3fffffffffffffff));
	std::vector<int64_t> result;
	result.reserve(count);
	for (size_t counter = 0; counter < count; ++counter) { result.push_back(int64_t(absDistr(absState))); }
	return result;
}

bool LongNumberTests::LongNumberAdditionTest::Run() const
{
	bool result = true;
	std::vector<int64_t> numbers = GetRandomNumbers(100);
	{ // both zeros
		result &= CheckPair(0, 0);
	}
	{ // one zero
		for each (auto number in numbers)
		{
			result &= CheckPair(0, number);
			result &= CheckPair(number, 0);
			result &= CheckPair(0, -number);
			result &= CheckPair(-number, 0);
		}
	}
	{ // non zero
		for each (auto left in numbers)
		{
			for each (auto right in numbers)
			{
				result &= CheckPair(left, right);
				result &= CheckPair(left, -right);
				result &= CheckPair(-left, right);
				result &= CheckPair(-left, -right);
			}
		}
	}
	{
		uint32_t nine = std::numeric_limits<uint32_t>::max();
		LongNumber ln1(std::vector<uint32_t>({ 1, nine, 345 }), false);
		LongNumber ln2(std::vector<uint32_t>({ 2, 0, 345 }), true);
		result &= (ln1 + ln2) == LongNumber(std::vector<uint32_t>({ 1, 0 }), true);
	}
	return result;
}

std::vector<int64_t> LongNumberTests::PerformanceCompare::GetRandomNumbers(const size_t & count) const
{
	std::random_device absSeed;
	std::mt19937_64 absState(absSeed());
	std::uniform_int_distribution<uint64_t> absDistr(uint64_t(0), uint64_t(0x3fffffffffffffff));
	std::random_device signSeed;
	std::mt19937_64 signState(signSeed());
	std::uniform_int_distribution<uint64_t> signDistr(uint64_t(0), uint64_t(1));

	std::vector<int64_t> result;
	result.reserve(count);
	for (size_t counter = 0; counter < count; ++counter)
	{
		result.push_back(int64_t(absDistr(absState)) * (int64_t(signDistr(signState)) * 2 - 1));
	}
	return result;
}

void LongNumberTests::PerformanceCompare::Run() const
{
	size_t N = 100;
	auto set0 = GetRandomNumbers(N);
	auto set1 = GetRandomNumbers(N);
	LongNumberIf lnif(int64_t(0));
	for (size_t i = 0; i < N; ++i)
	{
		LongNumberIf temp(set0[i]);
		temp += LongNumberIf(set1[i]);
		lnif += temp;
	}
	std::cout << lnif << std::endl;
	LongNumber ln(int64_t(0));
	for (size_t i = 0; i < N; ++i)
	{
		LongNumber temp(set0[i]);
		temp += LongNumber(set1[i]);
		ln += temp;
	}
	std::cout << ln << std::endl;
}

bool LongNumberTests::LongNumberMultiplicationTest::Run() const
{
	bool res = true;
	{
		uint32_t nine = std::numeric_limits<uint32_t>::max();
		LongNumber big(std::vector<uint32_t>({ nine, nine, nine, nine, nine, nine, nine }));
		LongNumber little(std::vector<uint32_t>({ nine, nine, nine, nine, nine }));
		LongNumber result(std::vector<uint32_t>({ nine, nine, nine, nine, nine - 1, nine, nine, 0, 0, 0, 0, 1 }));
		res &= (result == (big * little));
		res &= (result == (little * big));
	}
	{
		LongNumber ln(std::vector<uint32_t>({ 834568734, 934579234, 3495872934, 34957293, 0 }));
		LongNumber one(std::vector<uint32_t>({ 1 }));
		res &= (ln == (ln * one));
		res &= (ln == (one * ln));
		LongNumber zero(std::vector<uint32_t>({ 0 }));
		res &= (zero == (ln * zero));
		res &= (zero == (zero * ln));
	}
	{
		std::vector<LongNumber> partial;
		LongNumber hundredF(uint32_t(1));
		for (size_t i = 0; i < 10; ++i)
		{
			LongNumber acc(uint32_t(1));
			for (uint32_t j = 10 * i + 1; j <= 10 * (i + 1); ++j)
			{
				acc = acc * LongNumber(j);
				hundredF = hundredF * LongNumber(j);
			}
			partial.push_back(acc);
		}
		LongNumber M(uint32_t(1));
		for each (auto num in partial) { M = M * num; }
		res &= (M == hundredF);
	}
	{
		LongNumber pos(std::vector<uint32_t>({ 1, 1, 1 }), false);
		LongNumber neg(std::vector<uint32_t>({ 1, 1, 1 }), true);
		LongNumber respos(std::vector<uint32_t>({ 1, 2, 3, 2, 1 }), false);
		LongNumber resneg(std::vector<uint32_t>({ 1, 2, 3, 2, 1 }), true);
		res &= (respos == (pos * pos));
		res &= (respos == (neg * neg));
		res &= (resneg == (neg * pos));
		res &= (resneg == (pos * neg));
	}
	return res;
}

bool LongNumberTests::LongNumberDivisionTest::CheckPair(const LongNumber & first, const LongNumber & second) const
{
	bool result = true;
	LongNumber remainder;
	LongNumber quotient;
	try
	{
		quotient = first.QuotientRemainder(second, remainder);
		result &= ((second * quotient + remainder) == first);
		result &= (remainder < second.Absolute());
		result &= (LongNumber(uint32_t(0)) <= remainder);
	}
	catch (const LongNumber::DivideByZeroException& exc) { result &= (LongNumber(uint32_t(0)) == second); }
	try
	{
		quotient = second.QuotientRemainder(first, remainder);
		result &= ((first * quotient + remainder) == second);
		result &= (remainder < first.Absolute());
		result &= (LongNumber(uint32_t(0)) <= remainder);
	}
	catch (const LongNumber::DivideByZeroException& exc) { result &= (LongNumber(uint32_t(0)) == first); }
	return result;
}

bool LongNumberTests::LongNumberDivisionTest::CheckPairSigns(const LongNumber & first, const LongNumber & second) const
{
	bool result = true;
	result &= CheckPair(first, second);
	result &= CheckPair(first, -second);
	result &= CheckPair(-first, second);
	result &= CheckPair(-first, -second);
	return result;
}

LongNumber LongNumberTests::LongNumberDivisionTest::GetRandomNumber()
{
	std::uniform_int_distribution<uint32_t> absDistr(0, std::numeric_limits<uint32_t>::max());
	std::vector<uint32_t> result;
	size_t size = absDistr(absState) % 100;
	result.reserve(size);
	for (size_t counter = 0; counter < size; ++counter) { result.push_back(absDistr(absState)); }
	return LongNumber(result, false);
}

LongNumberTests::LongNumberDivisionTest::LongNumberDivisionTest()
{
	std::random_device absSeed;
	absState = std::mt19937_64(absSeed());
}

bool LongNumberTests::LongNumberDivisionTest::Run()
{
	uint32_t nine = std::numeric_limits<uint32_t>::max();
	bool res = true;
	LongNumber fb1002("113796925398360272257523782552224175572745930353730513145086634176691092536145985470146129334641866902783673042322088625863396052888690096969577173696370562180400527049497109023054114771394568040040412172632376");
	LongNumber fb501("225591516161936330872512695036072072046011324913758190588638866418474627738686883405015987052796968498626");
	res &= CheckPairSigns(fb1002, fb501);
	res &= CheckPairSigns(
		LongNumber("55325504782677307669912593071"), 
		LongNumber("12803182046945960413928791303578344629152701191943753805351675696880101068736172128056821462276339747823612636626337"));
	res &= CheckPairSigns(
		LongNumber(std::vector<uint32_t>({ nine, nine, nine, 0, 1, 1, nine })),
		LongNumber(std::vector<uint32_t>({ nine, nine, nine })));
	res &= CheckPairSigns(
		LongNumber(std::vector<uint32_t>({ nine, nine, nine, 0, 1, 1, nine })),
		LongNumber(std::vector<uint32_t>({ nine >> 7, nine, nine })));
	res &= CheckPairSigns(
		GetRandomNumber(), 
		LongNumber(uint32_t(1)));
	res &= CheckPairSigns(
		LongNumber(uint32_t(1)), 
		LongNumber(uint32_t(1)));
	res &= CheckPairSigns(
		LongNumber(uint32_t(0)), 
		LongNumber(uint32_t(0)));
	res &= CheckPairSigns(
		LongNumber(uint32_t(0)), 
		GetRandomNumber());
	for (int i = 0; i < 100; ++i) { res &= CheckPairSigns(GetRandomNumber(), GetRandomNumber()); }
	return res;
}

LongNumber LongNumberTests::LongNumberGCDTest::GCD(const LongNumber & first, const LongNumber & second) const
{
	// first >= 0 second >= 0
	bool f_less_s = first < second;
	LongNumber M = f_less_s ? second : first;
	LongNumber N = f_less_s ? first : second;
	LongNumber R;
	LongNumber* Mref = &M;
	LongNumber* Nref = &N;
	LongNumber* Rref = &R;
	while (*Nref != LongNumber(uint32_t(0)))
	{
		Mref->QuotientRemainder(*Nref, *Rref);
		std::swap(Mref, Nref);
		std::swap(Nref, Rref);
	}
	return *Mref;
}

bool LongNumberTests::LongNumberGCDTest::CheckTrio(const LongNumber & first, const LongNumber & second, const LongNumber & common) const
{
	bool result = true;
	LongNumber fsGCD = LongNumber::GCD(first, second);
	if (LongNumber(uint32_t(0)) != fsGCD)
	{
		LongNumber r1(uint32_t(1));
		LongNumber q1 = first.QuotientRemainder(fsGCD, r1);
		LongNumber r2(uint32_t(1));
		LongNumber q2 = second.QuotientRemainder(fsGCD, r2);
		result &= LongNumber(uint32_t(0)) == r1;
		result &= LongNumber(uint32_t(0)) == r2;
		result &= LongNumber(uint32_t(1)) == LongNumber::GCD(q1, q2);
	}
	//std::cout << "gcd" << fsGCD << std::endl;
	result &= GCD(first, second) == fsGCD;
	result &= LongNumber::GCD(second, first) == fsGCD;
	result &= (fsGCD * common) == LongNumber::GCD(first * common, second * common);
	return result;
}

LongNumber LongNumberTests::LongNumberGCDTest::GetRandomNumber()
{
	std::uniform_int_distribution<uint32_t> absDistr(0, std::numeric_limits<uint32_t>::max());
	std::vector<uint32_t> result;
	size_t size = absDistr(absState) % 10;
	result.reserve(size);
	for (size_t counter = 0; counter < size; ++counter) { result.push_back(absDistr(absState)); }
	return LongNumber(result, false);
}

LongNumberTests::LongNumberGCDTest::LongNumberGCDTest()
{
	std::random_device absSeed;
	absState = std::mt19937_64(absSeed());
}

bool LongNumberTests::LongNumberGCDTest::Run()
{
	bool res = true;
	uint32_t nine = std::numeric_limits<uint32_t>::max();
	LongNumber Nine(nine);
	LongNumber Zero(uint32_t(0));
	LongNumber One(uint32_t(1));
	LongNumber Even(uint32_t(34533444));
	LongNumber Ten(std::vector<uint32_t>({ 1, 0 }));
	res &= CheckTrio(GetRandomNumber() + Nine, GetRandomNumber() + One, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(GetRandomNumber() + Nine, GetRandomNumber() + One, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, One, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, One, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, Even, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, Even, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	res &= CheckTrio(Even, Ten, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(Even, Ten, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, Zero, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, Zero, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	res &= CheckTrio(Zero, Zero, LongNumber(std::vector<uint32_t>({ 1, 0, 0, 0, 0 }), false));
	res &= CheckTrio(Nine, Ten, LongNumber(std::vector<uint32_t>({ (34563 << 7), 0, 0, 0, 0 }), false));
	for (int i = 0; i < 100; ++i) { res &= CheckTrio(GetRandomNumber(), GetRandomNumber(), GetRandomNumber()); }
	return res;
}
