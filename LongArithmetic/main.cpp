#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include "LongNumber.h"
#include "LongNumberTests.h"


int main()
{
	auto start = std::chrono::steady_clock::now();
	LongNumber nine(uint32_t(9));
	LongNumber acc(uint32_t(1));
	for (size_t i = 0; i < 7000; ++i)
	{
		acc = acc * nine;
	}
	std::cout << acc << std::endl;
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start);
	std::cout << "Time: " << duration.count() << std::endl;

	LongNumber Na(uint32_t(1));
	LongNumber Nb(uint32_t(1));
	LongNumber Nc;
	for (size_t i = 0; i < 1000; ++i)
	{
		if (i == 500) { Nc = Na; }
		if (0 != (i & 1))
		{
			Nb += Na;
			//std::cout << Nb << std::endl;
		}
		else
		{
			Na += Nb;
			//std::cout << Na << std::endl;
		}
	}
	LongNumber fb1002("113796925398360272257523782552224175572745930353730513145086634176691092536145985470146129334641866902783673042322088625863396052888690096969577173696370562180400527049497109023054114771394568040040412172632376");
	LongNumber fb501("225591516161936330872512695036072072046011324913758190588638866418474627738686883405015987052796968498626");
	std::cout << "Nb: " << (Nb == fb1002) << std::endl;
	std::cout << "Nc: " << (Nc == fb501) << std::endl;

	LongNumberTests test;
	std::cout << "test: " << test.TestLongNumbers() << std::endl;

	int fg;
	std::cin >> fg;
	return 0;
}