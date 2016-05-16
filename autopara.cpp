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
#include <ctime>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstring>
#include<math.h>
#include <iostream>
#include <fstream>
#include <string>
# define TTLRange 10
# define beta 0.5
# define interval 1.5

using namespace std;
struct RulePrioSorter{
    bool operator ()(const array<int,2> &a, const array<int,2> &b) {
        return a[0]<b[0];
    }
};

void sortM(FlowRuleTable *oldMat)
{
    //   flowRuleTable newMat = *oldMat;
    int len = oldMat->get_rulenum();
    int t;
    for (int i=1;i<=len-1;i++)
    {
        for (int j=i+1;j<=len;j++)
        {
            if (oldMat->get(1,i) < oldMat->get(1,j))
            {
                t = oldMat->get(1,i);
                oldMat->set(1,i,oldMat->get(1,j));
                oldMat->set(1,j,t);
                t = oldMat->get(2,i);
                oldMat->set(2,i,oldMat->get(2,j));
                oldMat->set(2,j,t);
            }
        }
    }
    //return oldMat;
};

RulePrioSorter rulePrioSorter;
void Automatic::paraGenerate(int nFlow0, int nRule0, double alpha, float TTLMax, FlowRuleTable & table, floatCounter & flowPara, floatCounter & TTL, int & flowInterest){
    
    int flag = 0;
    int matchNum, maxRule, flag2;
    nFlow=nFlow0;
    nRule0=nRule;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_real_distribution<double> distDouble (0.0,1.0);
    uniform_int_distribution<int> distInt(1, TTLRange);
    FlowRuleTable oldSeq(2, nRule);
    // flowRuleTable newSeq(2, nRule);
    FlowRuleTable newPro(2, nRule);
    FlowRuleTable posSet(1,nFlow);
    int posLen=0;
    
    while (flag == 0)
    {
        table.setZero();
        for(int i=1; i<=nFlow; ++i)
        {
            for (int j=1; j<= nRule; ++j)
            {
                if (distDouble(generator) < beta)
                {
                    table.set(i,j,j);
                }
                
            }
            flowPara[i] = distDouble(generator);
        }
        flag = 1;
        for (int j = 1; j<=nRule; ++j)
        {
            if (table.col(j-1).sum() == 0)
            {
                flag = 0;
                continue;
            }
        }
        for (int j = 1; j<=nFlow; ++j)
        {
            if (table.row(j-1).sum() == 0)
            {
                flag = 0;
                continue;
            }
        }
        
        oldSeq.setZero();
        for (int i=1; i<= nRule; i++)
        {
            oldSeq.set(2,i,i);
            for (int j=1; j<=nFlow; j++)
            {
                if (table.get(j, i) > 0)
                {
                    oldSeq.set(1,i, 1+(oldSeq.get(1,i)));
                }
            }
        }
        
        sortM(&oldSeq);
        newPro.setZero();
        for (int i=1; i<=nRule; i++)
        {
            newPro.set(1,i,oldSeq.get(2,i));
            newPro.set(2,i,i);
        }
        sortM(&newPro);
        for (int i=1; i<=nRule; i++)
        {
            for (int j=1; j<=nRule; j++)
            {
                if (table.get(j,i) > 0)
                {
                    table.set(j,i,newPro.get(2,nRule-i+1));
                }
            }
        }
        
        posSet.setZero();
        posLen = 0;
        for (int j=1; j<=nFlow; j++)
        {
            matchNum = 0;
            maxRule = table.get_high_rule(j);
            for (int k=1; k<=nRule; k++)
            {
                if(table.get(j,k) > 0)
                {
                    matchNum++;
                }
            }
            if (matchNum >= 2)
            {
                flag2 = 0;
                for (int k = 1; k<=nFlow; k++)
                {
                    if ((k!=j) && (table.get(k,maxRule) > 0))
                    {
                        int mRule = 0;
                        for (int u = 1; u<= nRule; u++)
                        {
                            if ((table.get(k,u) > 0) && (table.get(j,u) == 0) && (table.get(k,u) < table.get(j,maxRule)))
                            {
                                mRule = 1;
                            }
                        }
                        if ((mRule == 0) && (poissonNumber(flowPara[k], 0, interval) > 0.5))
                        {
                            flag2 = 1;
                        }
                    }
                }
                if ((flag2 == 1) && (poissonNumber(flowPara[j], 0, interval) > 0.7) && (poissonNumber(flowPara[j],0,interval)<0.8))
                {
                    posLen++;
                    posSet.set(1,posLen,j);
                    
                    
                }
            }
        }
        
        if (posLen == 0)
        {
            flag = 0;
            // cout<<"len==0"<<endl;
            //cout<<table<<endl;
        }
        //cout<<"rubn"<<endl;
    }
    // cout<<"rubn"<<endl;
    int rd = ceil(distDouble(generator) * posLen);
    flowInterest = posSet.get(1,rd);
    
    for(int i=1;i<=nRule;++i){
        
        int rd = distInt(generator);
        TTL[i]=TTLMax * rd / TTLRange;
    }
    
    mSize = floor(alpha * nRule);
    if (mSize < 1)
    {
        mSize = 1;
    }
    // cout<<"rubn finish"<<endl;
    
}


void Automatic::save(string path){
    times++;
    ofstream myfile;
    myfile.open(path+"/para"+to_string(times)+".txt");
    myfile<<(*flowRuleTable);
    myfile<<endl;
    for(int i=0;i<nRule;i++)
        myfile<<TTL->get(i)<<"\t";
    myfile<<endl;
    for(int i=0;i<nFlow;i++)
        myfile<<flowPara->get(i)<<"\t";
    myfile<<endl<<"------target--------\n";
    myfile<<target<<endl;
    myfile<<"------result--------\n";
    for(std::set<int>::iterator it=attackFlow.begin(); it!=attackFlow.end(); ++it)
        myfile<<*it<<endl;
    myfile<<PrXQ.transpose()<<endl;
    myfile<<IG.transpose()<<endl;
    
    myfile.close();
}
int Automatic::generate(){
    //vector<> resultV = zeros(1,12);
    
    double TTLMax = 1;
    double alpha = 0.7;
    
    nFlow = 12;
    nRule = 10;
    int count = 120;
    int dCounter = 0;
    int i = 1;
    flowRuleTable->resize(nFlow, nRule);
    flowPara->resize(nFlow);
    TTL->resize(nRule);
    int flowInterest;
    qNum=1;
    double delta = 0.001;
    double limit = 0.001;
    //  double interval = 1.5;
    int maxm = 10000;
    for(int i=0;i<100;++i){
        paraGenerate(nFlow, nRule, alpha, TTLMax,*flowRuleTable,*flowPara, *TTL,flowInterest);
        cout<<*flowRuleTable<<endl;
        for(int i=0;i<nRule;i++)
            cout<<(*TTL)[i+1]<<" ";
        cout<<endl;
        for(int i=0;i<nFlow;i++){
            cout<<(*flowPara)[i+1]<<" ";
        }
        cout<<flowInterest<<endl;
        target=flowInterest;
        attackFlow.empty();
        updated=true;
        run(qNum,target,initialStateNum,attackFlow,PrXQ,IG);
        bool record_case=false;
        //save("/Users/ziqiaozhou/GoogleDrive/sdncode/database");
        if(attackFlow.count(target)==0){
            for(std::set<int>::iterator it=attackFlow.begin(); it!=attackFlow.end(); ++it){
                if (PrXQ((*it)-1,0)>0.5 && PrXQ((*it)-1,1)<1) {
                    record_case=true;
                    choose=*it;
                    break;
                }
            }
        }
        if(record_case)
            save("/Users/ziqiaozhou/GoogleDrive/sdncode/database");
        
    }
    
    return flowInterest;
}
