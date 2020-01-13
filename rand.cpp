#include <iostream>
#include <ctime>
#include <random>
using namespace std;

int unifRandInt(int a,int b) {
	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(a,b);
	return dis(this->gen);	
}
double unifRand() {
	return rand() / double (RAND_MAX);
}
double unifRand(double a, double b) {
	return (b-a)*unifRand() + a;
}

