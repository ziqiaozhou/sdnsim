//
//  attacker.h
//  sdnmodel
//
//  Created by ziqiao zhou on 5/13/16.
//  Copyright (c) 2016 ziqiao zhou. All rights reserved.
//

#ifndef sdnmodel_attacker_h
#define sdnmodel_attacker_h
#include "model3.h"
#include <vector>
#include<set>

class Attacker:public model3{
public:
    int qNum;
    int queryNum;
    Attacker(floatCounter * flowPara, FlowRuleTable *flowRuleTable, floatCounter * TTL, int mSize, int initialStateNum, double interval, double unit, double delta):model3(flowPara,flowRuleTable,TTL,mSize,initialStateNum,interval,unit,delta){
        //StateType i,j=1;
    };

    double flowProbCompute3M(StateProb2& stateProb0,std::set<int> &queryInterest,int queryResult);
    void conditionalEntropyComputeM3(int flowInterest,int initialStateNum,std::vector<double>&conditionalEntropyQ,std::vector<std::vector<double>>&PrXQ);
    void run(int qNum0,int flowInterest,StateType initialStateNum,std::set<int>&attackFlow, std::vector<std::vector<double>>*&PrXQ,std::vector<double>*&IG);
        
};
StateType num2StateAttack(int num, int flowNum,int qNum);
std::set<int> num2SetAttack(int num, int flowNum,int qNum);
#endif
