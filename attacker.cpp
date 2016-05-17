//
//  attacker.cpp
//  sdnmodel
//
//  Created by ziqiao zhou on 5/12/16.
//  Copyright (c) 2016 ziqiao zhou. All rights reserved.
//

#include "attacker.h"
#include <deque>
#include "sdnsim.h"
#include<cmath>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include<set>
#include<omp.h>
#include <unsupported/Eigen/MatrixFunctions>

#define entropy(pr) (((pr<=0)||(pr>=1))?0:(-pr*log(pr)-(1-pr)*log(1-pr)))
using namespace std;
//function [queryNum] = stateNumComputeAttack(flowNum, qNum)

//function [attackFlow, PrXQ, IG] = attackM3(flowPara, flowRuleTable, interval, flowInterest, initialStateNum, TTL, mSize, unit, delta, qNum)

StateType num2StateAttack(int num, int flowNum,int qNum){
    StateType state=0;
    int k = 1;
    int x,y;
    while (k <= qNum){
        x = factorial(flowNum - k) / factorial(flowNum - qNum);
        y = 0;
        while (num>x*y)
            ++y;
        state|=y?(1<<(y-1)):0;
        num = num - (y - 1) * x;
        k = k + 1;
    }
    return state;
}

set<int> bin2SetAttack(int num, int flowNum,int qNum){
    ++num;
    set<int> state;
    //StateType state=0;
    int k = 1;
    int x,y;
    while (k <= qNum){
        x = factorial(flowNum - k,flowNum - qNum);
        y = 0;
        while (num>x*y)
            ++y;
        state.insert(y);
        num = num - (y - 1) * x;
        k = k + 1;
    }
    return state;
}


double Attacker::flowProbCompute3M(StateProb2& stateProb0,set<int>& queryInterest,int queryResult){
    double queryProb = 0;
    StateType r = 0;
    //for(size_t i=0,row=stateProb0.rows();i<row;++i){
    for (StateProb2::InnerIterator it(stateProb0); it; ++it) {
        StateType oldList = legalState[it.index()];
        // int query=queryInterest;
        int t,count=-1;
        r = 0;
        for(std::set<int>::iterator it=queryInterest.begin(); it!=queryInterest.end(); ++it){
            t=*it;
            ++count;
            bool flag=false;
            for(int j=1;j<=nRule;++j){
                if (flowRuleTable->get(t, j) && exist_bit(oldList, j)){
                    r +=1<<count;
                    break;
                }
            }
        }
        if (r == queryResult){
            queryProb+= it.value();
        }
    }
    return queryProb;
}

void Attacker::conditionalEntropyComputeM3(int flowInterest,int initialStateNum,VecD &conditionalEntropyQ,MatD&PrXQ){
    int valueNum = pow(2,qNum);
    //  double ** PrQ=(double **)malloc(sizeof(double *)*queryNum);
    //for(int i=0;i<valueNum;++i){
    //  PrQ[i]=(double *)malloc(sizeof(double *)*valueNum);
    //}
    StateProb2 stateProbI(stateNum);
    stateProbI.reserve(stateNum);
    stateProbI.setZero();
    /* for (int i=0; i<stateNum; ++i) {
     stateProbI[legalState[i]]=0;
     }*/
    stateProbI.insert(initialStateNum) = 1;
    StateProb2 stateProbA=stateProbI;
    //stateProbA[initialStateNum] = 1;
    cout<<"trans start\n";
    TransProb TransA;
    /*if(Trans.size()==0)
     TransA =transComputation(flowInterest) ;
     else
     TransA = transComputation_ignore(flowInterest);*/
    /*if(updated){
     transComputation() ;
     updated=false;
     }
     cout<<"trans complete\n";
     TransA=transComputation_ignore(flowInterest);
     */
    TransA=transComputation(flowInterest);
   // cout<<Trans<<endl;
    //cout<<TransA<<endl;
    
    // cout<<Trans.nonZeros()<<endl;
    
    cout<<"trans A complete\n";
    cout<<TransA.nonZeros()<<endl;
    for (int i=0;i<fn;++i) {
        stateProbI=Trans*stateProbI;
        stateProbA=TransA*stateProbA;
    }
    cout<<"multiply\n";
    double PrQ;
    conditionalEntropyQ.setZero();
#pragma omp parallel for
    for (StateType i=0;i<queryNum;++i){
        set<int> state=bin2SetAttack(i,nFlow,qNum);
        for(int j=0;j<valueNum;++j){
            double PrQ=flowProbCompute3M(stateProbI, state, j);
            if (PrQ > 0){
                PrXQ(i,j)=flowProbCompute3M(stateProbA, state, j)/PrQ;
                conditionalEntropyQ(i)+=(PrQ*entropy(PrXQ(i,j)));
            }
        }
    }
    stateProb=stateProbI;
}

void Attacker::run(int qNum0,int flowInterest,StateType initialStateNum,set<int>&attackFlow, MatD &PrXQ,VecD &IG){
    init();
    updated=true;
    cout<<"run1";
    double pr = pow((1 - flowProb[flowInterest]), fn);
    /* for (int i=0;i<nFlow+1; ++i) {
     cout<<"flowprob="<<flowProb[i]<<endl;
     }*/
    double entropyQ = entropy(pr);
    double maxm = 0;
    qNum=qNum0;
    
    int flowNum=flowRuleTable->get_flownum();
    queryNum=((qNum > 0) && (flowNum >= qNum))?(factorial(flowNum,flowNum - qNum)):0;
    
    int valueNum = 1<<qNum;
    VecD conditionalEntropyQ(queryNum);//=(double *)zmalloc(sizeof(double)*queryNum);
    conditionalEntropyQ.setZero();
    PrXQ.resize(queryNum, valueNum);
    IG.resize(queryNum);
    IG.setOnes();
    conditionalEntropyComputeM3(flowInterest, initialStateNum,conditionalEntropyQ,PrXQ);
    IG=entropyQ*IG-conditionalEntropyQ;
    //IG.maxCoeff(&query,&tmp);
    for (int i = 0;i<queryNum;++i){
        if ( IG(i) > maxm  ){
            maxm = IG(i);
            attackFlow = bin2SetAttack(i, flowNum, qNum);
        }
        if(abs(IG(i) - maxm) <0.0000000001){
            set<int> attackFlowMid = bin2SetAttack(i, flowNum, qNum);
            if(attackFlowMid.count(flowInterest)){
                attackFlow = attackFlowMid;
                maxm = IG(i);
            }
        }
    }
    /*for(std::set<int>::iterator it=attackFlow.begin(); it!=attackFlow.end(); ++it){
     cout<<"choose"<<*it<<endl;
     }
     cout<<"choose end"<<endl;*/
}
