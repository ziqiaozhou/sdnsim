#include <iostream>
#include <list>
#include <numeric>
#include<iterator>

using namespace std;
typedef deque<int> LISTINT;
typedef deque<float> LISTINT;

LISTINT num2State3(int num, int ruleNum, int mSize);
unsigned nChoosek(unsigned n, unsigned k);
int nonZeroNum(LISTINT state);
LISTINT  num2bin(int stateNumber, int ruleNumber);


int main(void)
{
    LISTINT state;
   int it;
    state = num2State(5, 4, 3);
    for (it = 0; it < state.size(); it++)
    {
        //cout<<state[it]<<" ";
    }
    //cout<<endl;
    return 0;
}


LISTINT num2State3(int num, int ruleNum, int mSize)
{
    LISTINT cmp, state;
    int y, j, i, k, n_next, q, it;
    if (num > 1)
    else
    {
        j = 0;
        y = num - 1;
        while (y > 0)
        {
            num = y;
            j = j + 1;
            y = num - nChoosek(ruleNum, j);
        }
    }

    i = j;
    k = 0;
    while (num > k)
    {
        i = i + 1;
        cmp = num2bin(i, ruleNum);
        n_next = nonZeroNum(cmp);
        if (n_next == j)
        {
            k = k + 1;
        }
    }


    for (it = 0; it < cmp.size(); it++)
    {
        if (cmp[it] > 0)
        {
            state.push_back(it+1);
        }
    }


    if (state.size() > mSize)
    {
        //cout<<'larger than cache size!';
    }
    return state;

}

unsigned nChoosek(unsigned n, unsigned k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int nonZeroNum(LISTINT state)
{
    int ct = 0;
    int it;
    for (it = 0; it < state.size(); it++)
    {
        if (state[it] > 0)
        {
            ct++;
        }
    }
    return ct;
}



LISTINT  num2bin(int stateNumber, int ruleNumber)
{
    LISTINT state;
    int i;
    for (i = 1; i<= ruleNumber; i++)
    {
        if (stateNumber > 0)
        {
             state.push_front(stateNumber % 2);
             stateNumber = stateNumber/2;
        }
        else
        {
            state.push_front(0);
        }
    }
    return state;
}


