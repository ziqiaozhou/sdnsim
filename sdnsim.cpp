#include<list>
#include <algorithm>
#include <cmath>
#include <numeric>
#include<iterator>
#include "sdnsim.h"
#include<unordered_map>
#include<array>
#include<omp.h>
#include <Eigen/Sparse>
using namespace std;

struct SeqFlowSorter seqflowSorter;
void combineSeq(list<SeqFlow>& seqflow,list<double> seq2,int flow)
{
    for(LISTFLOAT::iterator it=seq2.begin();it!=seq2.end();it++)
    {
        seqflow.push_back(SeqFlow(*it,flow));
    }
    seqflow.sort(seqflowSorter);
    // sort(seqflow.begin(),seqflow.end(),seqflowSorter);
}
/* if (seq1(1,1)==-1)
 k=length(seq2);
 newSeq=zeros(k,2);
 for i=1:k
 newSeq(i,:)=[seq2(i),flow];
 end
 else
 newSeq=seq1;
 j=1;
 for i= 1: length(seq2)
 [k,~]=size(newSeq);
 while (( j <= k ) && ( newSeq(j,1)< seq2(i) ) )
 j=j+1;
 end
 if (j > k)
 newSeq(j,:)=[seq2(i),flow];
 j=j+1;
 else
 for t=k : -1 : j
 newSeq(t+1,:)=newSeq(t,:);
 end
 newSeq(j,:)=[seq2(i),flow];
 end
 
 end
 end
 
 end
 */

/*   newList = oldList;
 for i= length(newList) : -1 : 1
 if (newList(i) == d)
 newList(i)=[];
 end
 end
 end*/


/*
 function [num] = ceilM(TTL, unit, delta)
 % to compute ceil(TTL / unit) meanwhile addressing the error problem
 
 num = 0;
 sum = 0;
 flag = 0;
 
 while ((sum < TTL) && (flag == 0))
 num = num + 1;
 sum = sum + unit;
 if ((sum < TTL + delta) && (sum > TTL - delta))
 flag = 1;
 end
 end
 end
 */


/*
 function [bool] = doesMatch(flow, list, flowRuleTable)
 % find if there is a rule in list matching flow
 bool = 0;
 if (~isempty(list))
 for i = 1:length(list)
 if (flowRuleTable(flow, list(i)) ~= 0)
 if (bool == 0)
 bool = list(i);
 else
 if (flowRuleTable(flow, list(i)) > flowRuleTable(flow, bool))
 bool = list(i);
 end
 end
 end
 end
 end
 end
 */
int model::doesMatch(int flow,StateType stateNum)
{
    //  int rule=0;
    int tmp_prio;
    int prio=0;
    int ruleNo;
    int matched_rule=0;
    // //cout<<"doesMatch:cache state,flow"<<stateNum<<flow<<endl;
    if (stateNum==0) {
        return 0;
    }
    forAllRinS(ruleNo, stateNum){
            //  //cout<<"doesMatch:cache2"<<endl;
            tmp_prio=flowRuleTable->get(flow, ruleNo);
            //  //cout<<"doesMatch:cache3"<<endl;
            if(tmp_prio>prio)
            {
                //rule = ruleNo;
                prio=tmp_prio;
                matched_rule=ruleNo;
            }
        }
    return matched_rule;
};




int nonZeroNum(LISTINT state)
{
    int ct = 0;
    for (LISTINT::iterator it = state.begin(); it != state.end(); it++)
    {
        if (*it > 0)
        {
            ct++;
        }
    }
    return ct;
}


int  num2bin(StateType stateNumber, int ruleNumber,BoolState * state,int& nonZeroNum)
{
    nonZeroNum=0;
    int i;
    state->reset();
    stateNumber--;
    
    for (i = 1; stateNumber>0; ++i)
    {
        int tmp=stateNumber & 1;
        stateNumber=stateNumber>>1;
        state->set(i,tmp);
        //////cout<<stateNumber<<",tmp="<<tmp<<"i"<<i<<"state="<<state->get(i)<<endl;
        if(tmp)
            nonZeroNum++;
    }
    return nonZeroNum;
}


/*
 function [prob] = ruleEVT(rule, list, flowRuleTable, mSize, flowPara, TTL, unit, delta)
 
 [~, ruleNum] = size(flowRuleTable);
 if (length(list) > mSize)
 error = 'larger than cache size!'
 end
 
 if (~ifContain(rule, list))
 error = 'rule does not exist!'
 else
 if (length(list) < mSize)
 prob = 0;
 else
 prob = 0;
 for k = mSize: ceilM(TTL(rule), unit, delta)
 p = 1;
 for i = 1: ruleNum
 
 rPara = triggerFlowP(flowPara, i, flowRuleTable, list);
 if (~ifContain(i, list))
 p = p * poissonNumber0(rPara, 0, min([k * unit, TTL(i)]));
 else
 if (rule == i)
 p = p * poissonNumber0(rPara, 0, (k -1)* unit) * poissonNumber(rPara, 1, unit);
 else
 p = p * (1 - poissonNumber0(rPara, 0, min([TTL(i), k * unit])));
 end
 end
 end
 prob = prob + p;
 end
 
 
 end
 end
 
 end*/

double model::ruleEVT(int rule, StateType list,bool full)
{
    ////cout<<"ruleEVT"<<endl;
    int ruleNum=flowRuleTable->get_rulenum();
    double ttlrule=TTL->get(rule);
    double prob = 0,p;
    if (!exist_bit(list,rule ))
        cerr<< "rule does not exist!";
    else
    {
        if (!full)
            prob = 0;
        else
        {
            prob = 0;
            int maxv=ceilM(ttlrule, unit, delta);
            #pragma omp parallel for  reduction(+:prob)
            for (int k = mSize; k<=maxv; ++k)
            {
                p = 1;
                for(int i = 1; i<=ruleNum; ++i)
                {
                    bool exist=exist_bit(list,i);
                    double rPara = triggerFlowP( i, list,exist);
                    if (!exist)
                    {
                        p = p*poissonNumber0(rPara, 0, min(k*unit, TTL->get(i)));
                    }
                    else
                    {
                        if (rule == i){
                            p = p * poissonNumber0(rPara, 0, (k -1)* unit) * poissonNumber1(rPara, 1, unit);
                        }
                        else{
                            p = p * (1 - poissonNumber0(rPara, 0, min(TTL->get(i), k*unit)));
                        }
                        
                    }
                    
                }
                prob = prob + p;
            }
        }
    }
    ////cout<<"ruleEVT"<<endl;
    return prob;
}
/*
 function [lambda] = triggerFlowP(flowPara, rule, flowRuleTable, list)
 
 [flowNum, ruleNum] = size(flowRuleTable);
 fSet = [];
 k = 0;
 if (ifContain(rule, list))
 for i = 1:flowNum
 if (flowRuleTable(i, rule) ~= 0)
 flag = 0;
 for j = 1:length(list)
 if (flowRuleTable(i, list(j)) > flowRuleTable(i, rule))
 flag = 1;
 end
 end
 if (flag == 0)
 k = k + 1;
 fSet(k) = i;
 end
 end
 end
 else
 for i = 1: flowNum
 if (flowRuleTable(i, rule) ~= 0)
 flg = 0;
 for q = 1:length(list)
 if (flowRuleTable(i, list(q)) ~= 0)
 flg = 1;
 end
 end
 if (flg == 0)
 flag = 0;
 for j = 1:ruleNum
 if (flowRuleTable(i, j) > flowRuleTable(i, rule))
 flag = 1;
 end
 end
 if (flag == 0)
 k = k + 1;
 fSet(k) = i;
 
 end
 end
 end
 end
 end
 %fSet
 lambda = 0;
 for i = 1:length(fSet)
 lambda = lambda + flowPara(fSet(i));
 end
 
 end*/
double model::triggerFlowP(int rule,StateType state,bool exist)
{
    
    int flowNum=flowRuleTable->get_flownum();
    int ruleNum=flowRuleTable->get_rulenum();
     double lambda=0;
    int flag,flg;
   // LISTINT fSet;
    int it;
    ////cout<<rule<<exit<<endl;
    if (exist)
    {
        //////cout<<"triggerFlowP:yes"<<endl;
        for( int i = 1; i<=flowNum; ++i)
        {
            int rule_prio=flowRuleTable->get(i, rule);
            //////cout<<"triggerFlowP:yes:rule_prio"<<rule_prio<<endl;
            if (rule_prio != 0)
            {
                flag = 0;
                //////cout<<"triggerFlowP:yes:for"<<rule_prio<<endl;
                
                ////cout<<"triggerFlowP:yes"<<i<<","<<it<<endl;
                forAllRinS(it, state)
                {
                    //     //cout<<"triggerFlowP:yes"<<i<<","<<it<<endl;
                    //////cout<<"state *it="<<*it<<endl;
                    if(flowRuleTable->get(i, it) > rule_prio)
                    {
                        //////cout<<"triggerFlowP:yes:flag=1"<<rule_prio<<endl;
                        flag = 1;
                        break;
                    }
                    //   //cout<<"triggerFlowP:yes2"<<endl;
                    
                }
                
                if (flag == 0)
                {
                  //  fSet.push_back(i);
                    lambda +=flowPara->get(i);
                    ////cout<<"push back"<<i<<endl;
                }
            }
        }
        //////cout<<"triggerFlowP:yes:endfor"<<endl;
    }
    else
    {
        // //cout<<"triggerFlowP:no"<<endl;
        for (int i=1; i<=flowNum; ++i)
        {
            int rule_prio=flowRuleTable->get(i, rule);
            if (rule_prio != 0)
            {
                flg = 0;
                forAllRinS(it, state)
                {
                    ////cout<<i<<","<<it<<"flowrul="<<endl;
                    if (flowRuleTable->get(i, it) != 0)
                    {
                        flg = 1;
                        break;
                    }
                }
                if (flg == 0)
                {
                    flag = 0;
                    for(int j=1; j<=ruleNum; ++j)
                    {
                        
                        if (flowRuleTable->get(i, j) > rule_prio)
                        {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0)
                    {
                        lambda +=flowPara->get(i);
                       // fSet.push_back(i);
                        ////cout<<"else push back"<<i<<endl;
                    }
                }
            }
        }
    }
////cout<<state<<",lambda"<<lambda;
    return lambda;
}


/*void probNormalization(StateProb2& midStateProb)
 {
 double sum=0;
 // for(int i=0)
 for(size_t i=0;i<midStateProb.size();++i)
 {
 sum+=midStateProb[i];
 }
 if(sum==1)
 return;
 if(sum==0){
 cerr<<"==0"<<endl;
 return;
 }
 for(size_t i=0;i<midStateProb.size();++i)
 {
 midStateProb[i]/=sum;
 }
 }*/


double model::TTLProb(int rule, StateType state)
{
    int ruleNum=nRule;
    double prob,p;
    double total,lambda;
    int statesize=nonZeroNum(state);
    bool full=(statesize==mSize);
    double ttlrule=TTL->get(rule);
    if (interval < ttlrule)
    {
        prob = 0;
        return prob;
    }
    prob = 1;
    // double *lambdas=(double *)malloc(sizeof(double)*ruleNum);
    vector<double> lambdas(ruleNum);
    //memset(lambdas,1,sizeof(double)*ruleNum);
    //   //cout<<"lambdas"<<endl;
    // BoolState exists=(ruleNum);
    total = 1;
    if (full){
        total = 0;
        for(int i=1; i<=ruleNum; ++i)
        {
            bool existi=exist_bit(state,i);
            lambda =triggerFlowP(i,state,existi);// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            lambdas[i-1]=lambda;
            if (i == rule){
                prob *=(1 - poissonNumber0(lambda, 0, unit)) * poissonNumber0(lambda, 0, ttlrule - unit);
                //     cout<<"prob"<<prob<<endl;
            }
            else {
                if (existi){
                    prob *=(1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i))));
                    // cout<<"prob"<<prob<<endl;
                }
                else{
                    prob  *= poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                    //  cout<<"prob"<<prob<<endl;
                }
            }
        }
        //   cout<<"prob"<<prob<<endl;
        if(prob==0)
            return 0;
        int it;
        //  next_valid_rule(state>>=1,state);
        total=0;
        forAllRinS(it, state)
        {
            int ceil0=ceilM(min(ttlrule,TTL->get(it)), unit, delta);
            for (int k =statesize; k<=ceil0; ++k)
            {
                p = 1.0;
                for(int j=1; j<=ruleNum; ++j)
                {
                    lambda=lambdas[j-1];
                    if (!exist_bit(state,j))
                        p *= poissonNumber0(lambda, 0, min(k * unit, TTL->get(j)));
                    else
                    {
                        if (j == it)
                            p*=poissonNumber0(lambda, 0, (k - 1) * unit) * (1 - poissonNumber0(lambda, 0, unit));
                        else
                            p *=(1 - poissonNumber0(lambda, 0, min(k * unit, TTL->get(j))));
                    }
                    ////cout<<ceil0<<"k="<<"lambda"<<lambda<<"p="<<p<<endl;
                }
                
                total += p;
            }
        }
        //cout<<"lambdas3"<<endl;
    }else{
        total = 1;
        
        for(int i=1; i<=ruleNum; ++i)
        {
            bool existi=exist_bit(state,i);
            //////cout<<"ttlprob:for"<<endl;
            // //cout<<"lambdas2"<<i<<state<<existi<<endl;
            lambda =triggerFlowP(i,state,existi);// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            //  //cout<<"lambdas4"<<endl;
            lambdas[i-1]=lambda;
            /* if (lambda==0) {
             return 0;
             }*/
            //////cout<<"ttlprob:after triggerFlowP"<<endl;
            if (i == rule){
                prob *=(1 - poissonNumber0(lambda, 0, unit)) * poissonNumber0(lambda, 0, ttlrule - unit);
                if (existi)
                    total = total*(1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i))));
                else
                    total = total*poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                
            } else {
                double tmp=poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                if (existi){
                    prob *=tmp;
                    total *=tmp;
                }
                else{
                    prob *= tmp;
                    total *=tmp;
                }
            }
            //lambda = lambdas[i-1];//triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            ////cout<<"lamda="<<lambda<<"total="<<total<<"poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)))="<<poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)))
            //<<"1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)))"<<(1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i))))<<endl;
        }
        
    }
    ////cout<<"total="<<total<<endl;
    ////cout<<"prob="<<prob<<endl;
    //cout<<"prob="<<prob<<endl;
    if(prob==0)
        return 0;
    /*if (total==0){
        for(int k=0;k<ruleNum;++k)
            cout<<lambdas[k]<<" ";
        cerr<<"rule="<<ruleNum<<"full="<<full<<"total=="<<total<<"prob="<<prob<<endl;
        
        prob=0;
    }else*/
        prob /= total;
    // free(lambdas);
    return prob;
}

/*void model::matMultiply(TransProb& mat, StateProb2& oldprob)
 {
 ////////cout<<"size"<<mat.size()<<endl;
 
 StateProb2 newprob;//(oldprob.size(),0);
 double prob;
 #pragma omp parallel for
 for(StateType i=1;i<=stateNum;++i){
 StateType row=legalState[i];
 for (StateType j=1; j<=stateNum; ++j) {
 StateType col=legalState[j];
 prob=mat(row,col);
 //StateProb2::iterator old=oldprob.find(j);
 newprob[row]+=prob*oldprob[col];
 }
 
 }
 oldprob=newprob;
 }*/
/*function [time] = unitComputation(flowPara, delta, limit, TTL)
 % delta is the minimal time unit, while limit is the fidelity of normalize the probability.
 flowNum = length(flowPara);
 time = 0;
 sum = 1;
 while (sum > 1 - limit)
 time = time + delta;
 prob = 1;
 for i = 1:flowNum
 prob = prob * poissonNumber(flowPara(i), 0, time);
 end
 sum = prob;
 for i = 1:flowNum
 sum = sum + prob * poissonNumber(flowPara(i), 1, time) / poissonNumber(flowPara(i), 0, time);
 end
 
 end
 time = time + delta;
 flag = 0;
 while (flag == 0)
 time = time - delta;
 flag = 1;
 for i = 1:length(TTL)
 rm = TTL(i);
 while (rm > time)
 rm = rm - time;
 end
 if (rm <= time - limit)
 flag = 0;
 end
 end
 end
 end*/
double unitComputation(floatCounter *flowPara,double delta,double limit, floatCounter *TTL){
    int flowNum=flowPara->size();
    double time = 0;
    double sum = 1,prob,rm;
    while (sum > 1 - limit){
        time+=delta;
        prob=1;
        for (int i=1;i<=flowNum;++i){
            prob = prob * poissonNumber0(flowPara->get(i), 0, time);
        }
        sum = prob;
        for (int i=1;i<=flowNum;++i){
            sum = sum + prob * poissonNumber1(flowPara->get(i), 1, time) / poissonNumber0(flowPara->get(i), 0, time);
        }
    }
    time+=delta;
    int flag = 0;
    
    while(flag==0){
        time = time - delta;
        flag = 1;
        for (int i=1;i<=TTL->size();++i){
            rm = TTL->get(i);
            while (rm > time)
                rm = rm - time;
            if (rm <= time - limit)
                flag = 0;
        }
    }
    return time;
}
