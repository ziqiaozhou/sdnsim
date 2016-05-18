#include <iostream>
#include "sdnsim.h"
#include<list>
#include<deque>
#include"model3.h"
#include "attacker.h"
#include <set>
#include <vector>
#include<cmath>
#include<bitset>
#include "autopara.h"
/*StateType clear_bit(StateType u,int b){
 StateType t=(1<<(b-1));
 t=~t;
 return u&t;
 };
 StateType set_bit(StateType u,int b){
 StateType t=(1<<(b-1));
 return u|t;
 };*/

#define entropy(pr) ((pr<=0)||(pr>=1))?0:(-pr*log(pr)-(1-pr)*log(1-pr))
using namespace std;
int main(int argc, char *argv[])
{
    
    StateProb2 test(4);
    //cout<<test;
    long double pr=0.5;
    int flowNum=4,ruleNum=3;
    int mSize=1;
    long double msizerate=0.7;
    StateType initialStateNum=0;
    long double interval=2;
   	long double limit=0.001;
    long double delta=0.001;
    int target=1;
    int qnum=1;
    long double TTLmax=1;
    int runtimes=10;
    string label="";
    if(argc>1)
        label=argv[1];
    if(argc>2)
        flowNum=stoi(argv[2]);
    if(argc>3)
        ruleNum=stoi(argv[3]);
    if(argc>4)
        msizerate=stod(argv[4]);
    if(argc>5)
        TTLmax=stod(argv[5]);
    if(argc>6)
        interval=stod(argv[6]);
    if(argc>7)
        runtimes=stoi(argv[7]);
    FlowRuleTable table(flowNum,ruleNum);//flowNUm,ruleNUm
    floatCounter flowPara(flowNum);
    floatCounter TTL(ruleNum);
    
    long double unit=0;//=unitComputation(&flowPara, delta, limit, &TTL);
     Automatic a(label,&flowPara,&table,&TTL,mSize,initialStateNum,interval,unit,delta);
    
    int interestflow=a.generate(flowNum, ruleNum, msizerate, TTLmax,interval,runtimes);
    return 0;
}
