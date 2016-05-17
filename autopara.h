//
//  autopara.h
//  sdnmodel
//
//  Created by ziqiao zhou on 5/14/16.
//  Copyright (c) 2016 ziqiao zhou. All rights reserved.
//

#ifndef __sdnmodel__autopara__
#define __sdnmodel__autopara__

#include <stdio.h>
#include "attacker.h"
#include <array>
#include <set>
#include<string>
class Automatic:public Attacker{
public:
    int times;
    int target;
    int choose;
    int TTLmax;
    MatD PrXQ;
    std::string label;
    std::set<int>attackFlow;
    VecD IG;
    Automatic(std::string label0,floatCounter * flowPara, FlowRuleTable *flowRuleTable, floatCounter * TTL, int mSize, int initialStateNum, double interval, double unit, double delta):Attacker(flowPara,flowRuleTable,TTL,mSize,initialStateNum,interval,unit,delta){
         times=0;
         label=label0;
    };
    int paraGenerate(int flowNum, int ruleNum, double alpha, float TTLMax, FlowRuleTable & table, floatCounter & flowPara, floatCounter & TTL, int & flowInterest);
    int generate(int flowNum, int ruleNum, double alpha, float TTLMax,int interval,int runtimes);
    void save(std::string path);
};

typedef std::array<int, 2> RulePrio;
#endif /* defined(__sdnmodel__autopara__) */
