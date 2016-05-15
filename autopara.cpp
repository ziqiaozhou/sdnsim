//
//  autopara.cpp
//  sdnmodel
//
//  Created by ziqiao zhou on 5/14/16.
//  Copyright (c) 2016 ziqiao zhou. All rights reserved.
//

#include "autopara.h"
#include "sdnsim.h"
#include <cstdlib>
#include <ctime;
#include <vector>
#include <algorithm>
using namespace std;
struct RulePrioSorter{
    bool operator ()(const array<int,2> &a, const array<int,2> &b) {
        return a[0]<b[0];
    }
};
RulePrioSorter rulePrioSorter;
void Automatic::paraGenerate(int flowNum,int ruleNum,double alpha,double TTLMax){
    double TTLRange = 5;
    double belta = 0.3;
    int flag = 0;
    while (flag == 0){
       // flowRuleTable.resize(flowNum,ruleNum);
    retry:
        flag = 1;
        for(int i = 1;i<=nFlow;i++){
            for (int j=1; j<=nRule; j++) {
                srand(time(NULL));
                if (rand() < belta)
                    flowRuleTable->get(i,j)=j;
            }
            if (flowRuleTable->block(i,0,1,nRule).sum()==0) {
                flag=0;
                goto retry;
               // break;
            }
            flowPara->get(i) = rand();
        }
        for (int j=0;j<ruleNum;j++){
            if(flowRuleTable->block(0,j,nFlow,1).sum()==0){
                flag=0;
                goto retry;
            }
        }
        vector<array<int,2>> oldSeq(nRule);
        memset(oldSeq.data(),0,nRule*sizeof(int)*2);
        for (int i=1; i<=nRule; i++) {
            oldSeq[i][1]=i;
            for(int j = 1;j<=nFlow;j++){
                if (flowRuleTable->get(j, i) > 0)
                    oldSeq[i][0]+=1;
            }
        }
        vector<array<int,2>> newSeq=oldSeq;
        sort(newSeq.begin(),newSeq.end(),rulePrioSorter);
        for (int i=1; i<=nRule; i++) {
            newSeq[i][1] =i;
        }
         sort(newSeq.begin(),newSeq.end(),rulePrioSorter);
        for (int i=1; i<=nRule; i++) {
            for(int j = 1;j<=nFlow;j++){
                if (flowRuleTable->get(j,i) > 0)
                    flowRuleTable->get(j,i) = newSeq[ruleNum - i][1];
            }
        }
         for(int j = 1;j<=nFlow;j++){
            int  matchNum = 0;
            int maxRule = 0;
             for (int i=1; i<=nRule; i++){
                 
             }
         }
    }
    
}
void Automatic::automatic(){
    //vector<> resultV = zeros(1,12);
    
    double TTLMax = 1;
    double alpha = 0.7;
    
    int flowNum = 12;
    RuleType ruleNum = 10;
    int count = 120;
    int dCounter = 0;
    int i = 1;
    nFlow=flowNum;
    nRule=ruleNum;
    double delta = 0.001;
    double limit = 0.001;
    double interval = 1.5;
    int maxm = 10000;
    for(int i=0;i<100;i++){
        paraGenerate(flowNum, ruleNum, alpha, TTLMax);
    }
}