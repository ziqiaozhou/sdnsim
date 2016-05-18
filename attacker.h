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
    bool updated;
    Attacker(floatCounter * flowPara, FlowRuleTable *flowRuleTable, floatCounter * TTL, int mSize, int initialStateNum, long double interval, long double unit, long double delta):model3(flowPara,flowRuleTable,TTL,mSize,initialStateNum,interval,unit,delta){
        //StateType i,j=1;
        updated=true;
    };
    void run(int qNum0,int flowInterest,StateType initialStateNum,std::set<std::vector<int>>&attackFlow, MatD &PrXQ,VecD &IG);
     long double flowProbCompute3M(StateProb2& stateProb0,std::vector<int> queryInterest,int queryResult);
    void conditionalEntropyComputeM3(int flowInterest,int initialStateNum,VecD&conditionalEntropyQ,MatD&PrXQ);
    void future_conditionalEntropyComputeM3(int flowInterest,int initialStateNum,VecD &conditionalEntropyQ,MatD&PrXQ,StateProb2 &stateProbA,bool);
     
    // void run(int qNum0,int flowInterest,StateType initialStateNum,std::set<int>&attackFlow, MatD&PrXQ,VecD&IG);
    
};
StateType num2StateAttack(int num, int flowNum,int qNum);
std::vector<int> bin2SetAttack(int num, int flowNum,int qNum);
#endif
