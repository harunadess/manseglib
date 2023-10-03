#ifndef __TEST_UTIL_H__
#define __TEST_UTIL_H__

#include <iostream>

using namespace std;

void printBinary(double d)
{
    unsigned long l = *reinterpret_cast<unsigned long*>(&d);
    for(int i = 63; i >= 0; --i)
    {
        cout << ((l >> i) & 1);
        if(i == 63 || i == 52)
            cout << " ";
        else if(i == 32)
            cout << "|";
    }
    cout << endl;
}

// print binary of float
void printBinary(float f)
{
    unsigned int l = *reinterpret_cast<unsigned int*>(&f);
    for(int i = 31; i >= 0; --i)
    {
        cout << ((l >> i) & 1);
        if(i == 31 || i == 23)
            cout << " ";
    }
    cout << endl;
}

// print binary of int
void printBinary(int l)
{
    for(int i = 31; i >= 0; --i)
    {
        cout << ((l >> i) & 1);
        if(i % 8 == 0 && i > 0)
            cout << " ";
    }
	cout << endl;
}

#endif // __TEST_UTIL_H__