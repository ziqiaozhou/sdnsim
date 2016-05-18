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
int Automatic::paraGenerate(int nFlow0, int nRule0, long double alpha, float TTLMax, FlowRuleTable & table, floatCounter & flowPara, floatCounter & TTL, int & flowInterest){
    
    int flag = 0;
    int matchNum, maxRule, flag2;
    nFlow=nFlow0;
    nRule=nRule0;
    srand (time(NULL));
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_real_distribution<long double> distDouble (0.0,1.0);
    uniform_int_distribution<int> distInt(1, TTLRange);
    uniform_int_distribution<int> distIntFP(1, 1000);
    FlowRuleTable oldSeq(2, nRule);
    // flowRuleTable newSeq(2, nRule);
    FlowRuleTable newPro(2, nRule);
    FlowRuleTable posSet(1,nFlow);
    int posLen=0;
    
    while (flag == 0)
    {
    retry:
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
            flowPara[i] =(long double)(distIntFP(generator))/1000 ;
        }
        flag = 1;
        for (int j = 1; j<=nRule; ++j)
        {
            if (table.col(j-1).sum() == 0)
            {
                flag = 0;
                goto retry;
                // continue;
            }
        }
        for (int j = 1; j<=nFlow; ++j)
        {
            if (table.row(j-1).sum() == 0)
            {
                flag = 0;
                goto retry;
                //continue;
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
            for (int j=1; j<=nFlow; j++)
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
            goto retry;
            // //cout<<"len==0"<<endl;
            ////cout<<table<<endl;
        }
        ////cout<<"rubn"<<endl;
    }
    *flowRuleTable=table;
    // //cout<<"rubn"<<endl;
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
    return flowInterest;
    // //cout<<"rubn finish"<<endl;
    
}


void Automatic::save(string path){
    times++;
    //label2=chrono::system_clock::now().time_since_epoch().count();
    ofstream myfile;
    myfile.open(path+"/"+to_string(nFlow)+"_"+to_string(nRule)+"para"+label+"_"+to_string(times)+".txt");
    myfile<<nFlow<<"\t"<<nRule<<"\t"<<mSize<<"\t"<<TTLmax<<endl;
    myfile<<(*flowRuleTable);
    myfile<<endl;
    for(int i=0;i<nRule;i++)
        myfile<<TTL->get(i+1)<<"\t";
    myfile<<endl;
    for(int i=0;i<nFlow;i++)
        myfile<<flowPara->get(i+1)<<"\t";
    myfile<<endl<<"------target--------\n";
    myfile<<target<<endl;
    myfile<<"------result--------\n";
    for(std::set<set<int>>::iterator oneset=attackFlow.begin(); oneset!=attackFlow.end(); ++oneset){
        for(std::set<int>::iterator it=oneset->begin(); it!=oneset->end(); ++it){
            myfile<<*it<<"\t";
        }
        myfile<<endl;
    }
    myfile<<"------details--------\n";
    myfile<<PrXQ.transpose()<<endl;
    myfile<<IG.transpose()<<endl;
    
    myfile.close();
}


void Automatic::save_tmp(string path){
    tmp_times++;
    //label2=chrono::system_clock::now().time_since_epoch().count();
    ofstream myfile;
    myfile.open(path+"/"+to_string(nFlow)+"_"+to_string(nRule)+"para"+label+"_"+to_string(tmp_times)+".txt");
    myfile<<nFlow<<"\t"<<nRule<<"\t"<<unit<<"\t"<<mSize<<"\t"<<TTLmax<<endl;
    myfile<<(*flowRuleTable);
    myfile<<endl;
    for(int i=0;i<nRule;i++)
        myfile<<TTL->get(i+1)<<"\t";
    myfile<<endl;
    for(int i=0;i<nFlow;i++)
        myfile<<flowPara->get(i+1)<<"\t";
    myfile<<endl<<"------target--------\n";
    myfile<<target<<endl;
    myfile<<"------result--------\n";
    for(std::set<set<int>>::iterator oneset=attackFlow.begin(); oneset!=attackFlow.end(); ++oneset){
        for(std::set<int>::iterator it=oneset->begin(); it!=oneset->end(); ++it){
            myfile<<*it<<"\t";
        }
        myfile<<endl;
    }
    myfile<<"------details--------\n";
    myfile<<PrXQ.transpose()<<endl;
    myfile<<IG.transpose()<<endl;
    
    myfile.close();
}
void Automatic::save(string path,int& counter){
    counter++;
    //label2=chrono::system_clock::now().time_since_epoch().count();
    ofstream myfile;
    myfile.open(path+"/"+to_string(nFlow)+"_"+to_string(nRule)+"para"+label+"_"+to_string(counter)+".txt");
    myfile<<nFlow<<"\t"<<nRule<<"\t"<<unit<<"\t"<<mSize<<"\t"<<TTLmax<<endl;
    myfile<<(*flowRuleTable);
    myfile<<endl;
    for(int i=0;i<nRule;i++)
        myfile<<TTL->get(i+1)<<"\t";
    myfile<<endl;
    for(int i=0;i<nFlow;i++)
        myfile<<flowPara->get(i+1)<<"\t";
    myfile<<endl<<"------target--------\n";
    myfile<<target<<endl;
    myfile<<"------result--------\n";
    for(std::set<set<int>>::iterator oneset=attackFlow.begin(); oneset!=attackFlow.end(); ++oneset){
        for(std::set<int>::iterator it=oneset->begin(); it!=oneset->end(); ++it){
            myfile<<*it<<"\t";
        }
        myfile<<endl;
    }
    myfile<<"------details--------\n";
    myfile<<PrXQ.transpose()<<endl;
    myfile<<IG.transpose()<<endl;
    myfile.close();
}

int Automatic::generate(int flowNum, int ruleNum, long double alpha, float TTLMax0,double interval0,int runTimes){
    //vector<> resultV = zeros(1,12);
    interval=interval0;
    nFlow = flowNum;
    nRule = ruleNum;
    TTLmax=TTLMax0;
    flowRuleTable->resize(nFlow, nRule);
    flowRuleTable->setZero();
    flowPara->resize(nFlow);
    flowPara->reset();
    TTL->resize(nRule);
    TTL->reset();
    int flowInterest;
    qNum=1;
    int counter=0;
    delta = 0.001;
    limit = 0.001;
    int i=0;
    //  long double interval = 1.5;
    // int maxm = 10000;
    while (i<runTimes){
        int tmp=paraGenerate(nFlow, nRule, alpha, TTLmax,*flowRuleTable,*flowPara, *TTL,flowInterest);
        ////cout<<"generated"<<endl;
        /*   //cout<<*flowRuleTable<<endl;
         for(int i=0;i<nRule;i++)
         //cout<<(*TTL)[i+1]<<" ";
         //cout<<endl;
         for(int i=0;i<nFlow;i++){
         //cout<<(*flowPara)[i+1]<<" ";
         }
         //cout<<flowInterest<<endl;*/
        target=flowInterest;
        attackFlow.empty();
        updated=true;
        
        run(qNum,target,initialStateNum,attackFlow,PrXQ,IG);
       // //cout<<"FP:"<<flowProb<<endl;
       // //cout<<"Trans:"<<Trans<<endl;
        bool record_case=false;
        //save("/Users/ziqiaozhou/GoogleDrive/sdncode/database");
        switch (caseType) {
            case CASE_SPECIAL:
                save("../special/",counter);
                break;
            default:
            {
                
                for(std::set<set<int>>::iterator oneset=attackFlow.begin(); oneset!=attackFlow.end(); ++oneset){
                    if(oneset->count(target)>0){
                        record_case=false;
                        break;
                    }
                    for(std::set<int>::iterator it=oneset->begin(); it!=oneset->end(); ++it){
                        if (PrXQ((*it)-1,0)>0.5 && PrXQ((*it)-1,1)<0.5) {
                            record_case=true;
                            choose=*it;
                            break;
                        }
                    }
                }
                if(record_case)
                   { save("../data/");
                    i = i + 1; }
                else
                {
                    save_tmp("../tmp/");
                    i = i + 1;}
            }
                break;
        }
    }
    return 0;
    //return flowInterest;
}
