#include"sdnsim.h"
#include<cmath>
#include<list>
#include <numeric>
#include"model3.h"
#include<omp.h>
#include <deque>
#include<set>

#include <unsupported/Eigen/MatrixFunctions>
using namespace std;

//deque<double> stateProb;

/*
 void model3::num2State(StateType num, int ruleNum, int mSize,LISTINT &state,BoolState& cmp)
 {
 
 StateType y, i, k;
 int nonZeroNum, q,j;
 if (num == 1)
 return;
 // //////cout<<"noe return"<<endl;
 j = 0;
 y = num - 1;
 
 
 while (y > 0&&j<ruleNum)
 {
 num = y;
 j = j + 1;
 y = num - nChoosek(ruleNum, j);
 //        //////cout<<num<<","<<"j="<<j<<"y="<<y<<endl;
 }
 i = j;
 k = 0;
 // BoolState cmp(ruleNum);
 while (num > k)
 {
 cmp.reset();
 i = i + 1;
 num2bin(i, ruleNum,&cmp,nonZeroNum);
 
 if (nonZeroNum == j)
 {
 k = k + 1;
 }
 }
 //	allStateMap[num]=cmp;
 
 for (int it = 1; it < cmp.size()+1; it++)
 {
 if (cmp[it] > 0)
 {
 //////cout<<"state add"<<it<<endl;
 state.push_back(it);
 }
 }
 // //////cout<<"statesize="<<state.size()<<endl;
 
 }*/

StateType model3::bin2Num(StateType state)
{
    ////cout<<"in state2Num3"<<endl;
    // size_t len=state.size();
    /*if (len > mSize)
     {
     cerr<< "rule list larger than cache size!"<<endl;
     return -1;
     }*/
    int num = 0;
    
    if (state > 0)
    {
        int len=nonZeroNum(state);
        StateType i = 1;
        for (i=1;i < len;i++){
            num = num + nChoosek(nRule, i);
            //i = i + 1;
        }
        // num=stateNumCompute3(ruleNum,mSize);
        int k = 0;
        int j = 1;
        while (state&& ((nRule-j)>=( len - k)) && (j < nRule))
        {
            if (state&1)
            {
                num+= nChoosek(nRule - j, len - k);
                k = k + 1;
            }
            state>>=1;
            j = j + 1;
        }
        num++;
    }
    ////cout<<"end state2Num3"<<endl;
    return num;
};


TransProb model3::transComputation(int ignored_flow)
{
    //omp_set_num_threads(2);
    
    //   TransProb Trans;
    Trans.resize(stateNum,stateNum);
    Trans.reserve(stateNum*stateNum/2);
    Trans.setZero();
    TransProb TransA=Trans;
    //#pragma omp parallel for
   // allnewStateProb.resize()
    for( StateType k=0;k<stateNum;k++){
      
        //   #pragma omp parallel for
        for (int j=0; j<=nFlow; j++)
        {
            StateProb2 newStateProb(stateNum,1);
            // newStateProb.setZero();
            //  //cout<<"newStateProb"<<endl;
            flowState(j, k,newStateProb);
            newStateProb*=flowProb[j];
         
            //#pragma omp critical
            {
                Trans.col(k)+=newStateProb;
                if(ignored_flow!=j)
                    TransA.col(k)+=newStateProb;
            }
        
        }
    }
    Trans.prune(0,0);
    TransA.prune(0,0);
    return TransA;
}



TransProb model3::transComputation()
{
    //omp_set_num_threads(2);
    int flowNum=nFlow;
    //   TransProb Trans;
    Trans.resize(stateNum,stateNum);
    Trans.reserve(stateNum*stateNum/2);
    Trans.setZero();
    //  TransProb::iterator it_in_trans;
    //StateType all=1<<ruleNum;
    //#pragma omp parallel for
    allnewStateProb.resize(nFlow);
  //  allnewStateProb.resize(stateNum,stateNum*nFlow);
    //allnewStateProb.reserve(stateNum*nRule*nFlow);
    for( StateType k=0;k<stateNum;k++){
        // StateType i=legalState[k];
     //   #pragma omp parallel for
        for (int j=0; j<=nFlow; j++)
        {
            if(allnewStateProb[j].size()<stateNum)
                allnewStateProb[j].resize(stateNum);
            StateProb2 newStateProb(stateNum,1);
           // newStateProb.setZero();
            //  //cout<<"newStateProb"<<endl;
            flowState(j, k,newStateProb);
         //  if(allnewStateProb.size()>j)
          //      allnewStateProb.block(0, j+k*nFlow, stateNum, 1)=newStateProb;
                allnewStateProb[j][k]=(newStateProb);
            
            // //cout<<"for"<<newStateProb.size()<<endl;
            //for( StateType k=0;k<legalState.size();k++){
           // cout<<"flow="<<j<<endl;
         //   cout<<"oldstate"<<legalState[k]<<endl;
           // cout<<newStateProb<<endl;
            newStateProb*=flowProb[j];
          //  cout<<"newStateProb"<<newStateProb<<endl;
           // block<stateNum,1>(1,k)
//#pragma omp critical
            {
            Trans.col(k)+=newStateProb;
            }
           /* for(size_t row = 0, nRows = newStateProb.rows(); row < nRows; ++row)
            {
                //cout<<i<<j<<"newStateProb"<<stateprob.first<<","<<stateprob.second<<endl;
                //if(newStateProb.find(k))
                //stateprob.second*=flowProb[j];
                // newStateProb[k]*=flowProb[j];
                //#pragma critical
                {
                    //      //////cout<<"check"<<"="<<stateprob.first<<","<<stateprob.second<<"i="<<i<<endl;
                    // it_in_trans=Trans.InTrans(stateprob.first,i);
                    Trans(row,k)+=newStateProb(row);
                    ////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
                }
            }*/
        }
    }
    Trans.prune(0,0);
    return Trans;
}

TransProb model3::transComputation_ignore(int ignored_flow)
{
    //omp_set_num_threads(2);
   /* if(Trans.size() ==0){
        transComputation();
    }*/
    TransProb newtrans=Trans;
    //TransProb::iterator it_in_trans;
#pragma omp parallel for
    for( StateType k=0;k<legalState.size();k++){
        // StateType i=legalState[k];
        int j=ignored_flow;
        StateProb2 newStateProb=allnewStateProb[j][k];
        double baseprob=flowProb(j);
        //  //cout<<"newStateProb"<<endl;
        
        // //cout<<"for"<<newStateProb.size()<<endl;
        //for( StateType k=0;k<legalState.size();k++){
        newStateProb*=baseprob;
        // newtrans=(newtrans-newStateProb);
        newtrans.col(k)-=newStateProb;
      /*  for(size_t row = 0, nRows = newStateProb.rows(); row < nRows; ++row)
        {
            //cout<<i<<j<<"newStateProb"<<stateprob.first<<","<<stateprob.second<<endl;
            //if(newStateProb.find(k))
            //  stateprob.second*=baseprob;
            // newStateProb[k]*=flowProb[j];
#pragma critical
            {
                //      //////cout<<"check"<<"="<<stateprob.first<<","<<stateprob.second<<"i="<<i<<endl;
                
                //it_in_trans=newtrans.InTrans(stateprob.first,i);
                //double now=newtrans[it_in_trans->first];
                ////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
                if (newtrans(row,k))
                {
                    ////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
                    newtrans(row,k)-=newStateProb(row);
                    //it_in_trans->second -=stateprob(i);
                    //  //////cout<<"it_in_trans"<<"="<<it_in_trans->second<<endl;
                }
            }
        }*/
    }
    return newtrans;
}


void model3::flowState(int flow, StateType oldStateNum,StateProb2& newStateProb)
{
        //    //////cout<<"flowState3:enter"<<endl;
    StateType oldbinState=legalState[oldStateNum];
   // cout<<flow<<"<-flow"<<oldbinState<<"<-oldbinary"<<endl;

    newStateProb.resize(stateNum, 1);
    newStateProb.setZero();
    int oldsize=nonZeroNum(oldbinState);
    newStateProb.reserve(((oldsize>= mSize)&&(flow > 0))?oldsize:1);
    int ruleNum=flowRuleTable->get_rulenum();
    LISTINT oldList,midList;
    // BoolState bool_oldList(ruleNum);
    double p;
    //  num2State3(oldStateNum, ruleNum, mSize,oldList,bool_oldList);
    //    //////cout<<"flowState3:num2State3"<<endl;
        int matched_prio=0;//matched_rule=0;
    int matched_rule=0;
    int high_prio_rule=flowRuleTable->get_high_rule(flow);
    // LISTINT::iterator matched_pos=oldList.end();
   // cout<<"stateNum="<<stateNum<<endl;
    StateProb2 midStateProb(stateNum),ttlStateProb(stateNum);
    midStateProb.setZero();
    ttlStateProb.setZero();
    size_t ttl_firststate=-1;
    // newStateProb.clear();
    // LISTINT newList;
    //  BoolState bool_newlist(ruleNum);//check whether a rule is in newlist in O(1)
    //  //cout<<"flowState3:cache"<<endl;
    
    //  //cout<<"flowState3:cache1"<<endl;
    matched_rule=doesMatch(flow,oldbinState);
    if ((oldsize>= mSize) && (flow > 0) )
    {
        
        if(!matched_rule){
            ////cout<<"flowState3:if"<<endl;
            //LISTINT::iterator it=oldList.begin();
            int deletedrule;
            StateType total=oldbinState;
            forAllRinS(deletedrule,total){
                StateType oneBinNum=clear_bit(oldbinState,deletedrule);
                StateType newBinNum=cacheLRU(flow, oneBinNum,(matched_rule!=deletedrule)?matched_rule:0);
                midStateProb.insert(bin2Num(newBinNum))=ruleEVT(deletedrule, oldbinState,true);
            }
            midStateProb/=midStateProb.sum();
        }else{
            StateType newBinNum=cacheLRU(flow, oldbinState ,matched_rule);
            midStateProb.setZero();
            midStateProb.insert(bin2Num(newBinNum))=1;
        }
        // //////cout<<"flowState3:norm"<<endl;
    }
    else
    {
        StateType newBinNum;
        newBinNum=cacheLRU(flow, oldbinState,matched_rule);
        //    //cout<<"flowState3:cache2"<<endl;
        midStateProb.setZero();
       midStateProb.insert(bin2Num(newBinNum))=1;
    }
    for (StateProb2::InnerIterator it(midStateProb); it; ++it) {
        StateType newList=legalState[it.index()];
        int k;
        int newrule=0;
        
        if (flow == 0)
            k = 1;
        else{
            k = 2;
            if (exist_bit(newList,high_prio_rule)) {
                newrule=high_prio_rule;
            }else{
                newrule=matched_rule;
            }
        }
        ttlStateProb.setZero();
        if (nolessbit(newList,k ))
        {
            //  //cout<<"noless"<<endl;
            
            ttlStateProb.insert(it.index())=it.value();
            ttl_firststate=it.index();
            int j_it;
            forAllRinS(j_it, newList)
            {
                if(j_it==newrule){
                    continue;
                }
               // cout<<"ttl table"<<j_it<<"index="<<it.index()<<ttlProbTable;
                p=ttlProbTable.coeffRef(it.index(),j_it-1);
                if(p){
                    p=(p>1)?0:p;
                }else{
                    p = TTLProb(j_it, newList);
                    ttlProbTable.insert(it.index(),j_it-1)=p?p:1.1;
                }
                if (p > 0)
                {
                    StateType onelist=clear_bit(newList, j_it);
                    double pp=p*it.value();
                    ttlStateProb.insert(bin2Num(onelist))=pp;
                    ttlStateProb.coeffRef(ttl_firststate)-=pp;
                }
            }
            //     //cout<<"flowState3:cache4"<<endl;
        }else
        {
            ttlStateProb.coeffRef(it.index())=it.value();
        }
        newStateProb+=ttlStateProb;
    }
    newStateProb.prune(0);
    //cout<<"newstateprob"<<newStateProb<<endl;
}

void model::initFlowProb(){
    
    if(!flowProb.size()){
       flowProb.resize(stateNum);
    }
    flowProb.resize(stateNum);
   
    double prob=1;
#pragma omp parallel for
    for (int i = 1; i<=nFlow; i++)
    {
        prob *= poissonNumber(flowPara->get(i), (int)0, unit);
    }
    
    //////cout<<"hello"<<endl;
    flowProb(0)=prob;
#pragma omp parallel for
    for(int i = 1; i<=nFlow; i++)
    {
        std::cout<<i<<std::endl;
        std::cout<<flowPara->get(i);
        flowProb(i) = poissonNumber(flowPara->get(i), (int)1, unit) * prob/ poissonNumber(flowPara->get(i), (int)0, unit);
    }
    flowProb/=flowProb.sum();
    //probNormalization(flowProb,1+nFlow);
}


void model3::run()
{
    //////cout<<"hello"<<endl;
    // int flowNum=flowRuleTable->get_flownum();
    //int ruleNum=flowRuleTable->get_rulenum();
    //////cout<<"hello"<<endl;
    // StateType stateNum = stateNumCompute(mSize,nRule);
    //////cout<<stateNum<<endl;
    //////cout<<"hello"<<endl;
    // int fn = ceilM(interval, unit, delta);
    
    // stateProb;
    // return stateProb;
    stateProb.resize(stateNum,1);
    stateProb.reserve(stateNum);
    stateProb.setZero();
    stateProb.insert(initialStateNum) = 1;
    
    //deque<double> flowProb(1+nFlow,1);
    // = (double*)malloc((1+flowNum)*sizeof(double));
    
    
    cout<<"build flowProb finished"<<endl;
    //////cout<<"hello"<<endl;
    TransProb Trans = transComputation();
    //////cout<<"hello end trans"<<endl;
    // double sumtrans[Trans.size()];
    // memset(sumtrans)
    for(size_t i=0,nrow=Trans.rows();i<nrow;i++){
        for(size_t j=0,ncols=Trans.cols();j<ncols;j++){
        //cout<<"trans";
        cout<<Trans.coeffRef(i,j)<<" ";
        }
        cout<<endl;
        //   sumtrans[ontran.first[1]-1]+=ontran.second;
    }
    cout<<"build trans finished"<<endl;
    for (int i=0; i<fn; i++) {
         stateProb=Trans*stateProb;
    }
   
    /*for(int i = 1;i<=fn;i++){
        //         //////cout<<"matMultiply"<<endl;
        matMultiply(Trans, stateProb);
    }*/
    //  //////cout<<"matMultiply end"<<endl;
    //////cout<<"hello statenhum="<<stateNum<<endl;
    //return stateProb;
}
