#include"sdnsim.h"
#include<cmath>
#include<list>
#include <numeric>
#include"model3.h"
#include<omp.h>
#include <deque>
#include<set>
#include <chrono>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace std::chrono;
//deque<long double> stateProb;
long double model3::TTLStateProb( StateType state,vector<long double> lambdas)
{
    long double maxt = 0,tmp;
    int ruleNo;
    forAllRinS(ruleNo,state){
        tmp=TTL->get(ruleNo);
        if(tmp>maxt){
            maxt=tmp;
        }
    }
    long double prob;
    long double total;
    int statesize=nonZeroNum(state);
    bool full=(statesize==mSize);
    prob = 1;
    
    if (full){
        total = 0;
        
        int ceil0=ceilM(maxt, unit, delta);
//#pragma omp parallel for reduction(+:total)
        for (int k =statesize; k<=ceil0; ++k)
        {
            int it;
            forAllRinS(it, state)
            {
                if ((TTL->get(it)> k * unit)){
                    long double lambda =lambdas[it-1];
                    long double p = poissonNumber0(lambda, 0, (k - 1) * unit) * (1 - poissonNumber0(lambda, 0, unit));
                    for(int j=1; j<=nRule; ++j)
                    {
                        if(j!=it){
                            lambda=lambdas[j-1];
                            if (!exist_bit(state, j)){
                                p *= poissonNumber0(lambda, 0, k * unit);
                            }
                            else{
                                p*= (1 - poissonNumber0(lambda, 0,  min(k * unit,TTL->get(j) ) ) );
                                
                            }
                        }
                        
                    }
                    
                    total+=p;
                    
                }
                
            }
        }
    }else{
        total=1;
//#pragma omp parallel for reduction(*:total)
        for(int i=1; i<=nRule; ++i)
        {
            long double lambda =lambdas[i-1];// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            if (exist_bit(state, i)){
                total *=(1 - poissonNumber0(lambda, 0, TTL->get(i)));
                // //cout<<"prob"<<prob<<endl;
            }
            else{
                total*= poissonNumber0(lambda, 0, TTL->get(i));
                //  //cout<<"prob"<<prob<<endl;
            }
        }
        
    }
    return total;
}

/*
 void model3::num2State(StateType num, int ruleNum, int mSize,LISTINT &state,BoolState& cmp)
 {
 
 StateType y, i, k;
 int nonZeroNum, q,j;
 if (num == 1)
 return;
 // ////////cout<<"noe return"<<endl;
 j = 0;
 y = num - 1;
 
 
 while (y > 0&&j<ruleNum)
 {
 num = y;
 j = j + 1;
 y = num - nChoosek(ruleNum, j);
 //        ////////cout<<num<<","<<"j="<<j<<"y="<<y<<endl;
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
 ////////cout<<"state add"<<it<<endl;
 state.push_back(it);
 }
 }
 // ////////cout<<"statesize="<<state.size()<<endl;
 
 }*/

StateType model3::bin2Num(StateType state)
{
    //////cout<<"in state2Num3"<<endl;
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
        for (i=1;i < len;++i){
            num += nChoosek(nRule, i);
            //i = i + 1;
        }
        // num=stateNumCompute3(ruleNum,mSize);
        int k = 0;
        int j = 1;
        while (state && (j < nRule)&& ((nRule-j)>=(len-k)))
        {
            if (state&1)
            {
                num+= nChoosek(nRule - j, len - k);
                ++k;
            }
            state>>=1;
            ++j;
        }
        num++;
    }
    //////cout<<"end state2Num3"<<endl;
    return num;
};


void model3::transComputation(int ignored_flow,TransProb & TransA)
{
    //omp_set_num_threads(2);
    
    //   TransProb Trans;
    Trans.resize(stateNum,stateNum);
    Trans.reserve(maxTrans);
   // //cout<<"hi";
    Trans.setZero();
    ////cout<<"hi2";
     TransA=Trans;
    
    //#pragma omp parallel for
    StateProb2 newStateProb(stateNum);
    for( StateType k=0;k<stateNum;++k){
        StateProb2 deletenewStateProb;
        
        //#pragma omp parallel for
        for (int j=0; j<=nFlow; ++j)
        {
            
            //newStateProb.reserve();
            //newStateProb.setZero();
            //  ////cout<<"newStateProb"<<endl;
           // //cout<<k;
            flowState(j, k,newStateProb);
            newStateProb*=flowProb[j];
            if(j==ignored_flow)
                deletenewStateProb=newStateProb;
            //#pragma omp critical(section1)
            {
                Trans.col(k)+=newStateProb;
            }
        }
        //#pragma omp critical(section1)
        {
            TransA.col(k)=Trans.col(k)-deletenewStateProb;
        }
        
    }
    Trans.prune(0,0);
    TransA.prune(0,0);
    //return TransA;
}



void model3::transComputation()
{
    //omp_set_num_threads(2);
    int flowNum=nFlow;
    //   TransProb Trans;
    Trans.resize(stateNum,stateNum);
    Trans.reserve(stateNum*nChoosek(nRule, nRule/2));
    Trans.setZero();
    //  TransProb::iterator it_in_trans;
    //StateType all=1<<ruleNum;
    //#pragma omp parallel for
    //  allnewStateProb.resize(stateNum,stateNum*nFlow);
    allnewStateProb.resize(stateNum);
    //allnewStateProb.reserve(stateNum*nRule*nFlow);
    StateProb2 newStateProb(stateNum,1);
    for( StateType k=0;k<stateNum;++k){
        // StateType i=legalState[k];
        //   #pragma omp parallel for
        
        allnewStateProb[k].resize(nFlow);
        for (int j=0; j<=nFlow; ++j)
        {
            allnewStateProb[k][j].resize(stateNum, 1);
            allnewStateProb[k][j].setZero();
            allnewStateProb[k][j].reserve(mSize);
            //  //cout<<k<<endl;
            // newStateProb.setZero();
            //  ////cout<<"newStateProb"<<endl;
            flowState(j, k,allnewStateProb[k][j]);
            
            allnewStateProb[k][j]*=flowProb[j];
            
            //  //cout<<"newStateProb"<<newStateProb<<endl;
            // block<stateNum,1>(1,k)
            //#pragma omp critical
            {
                Trans.col(k)+=allnewStateProb[k][j];
            }
            // allnewStateProb[k][j]=(newStateProb);
            /* for(size_t row = 0, nRows = newStateProb.rows(); row < nRows; ++row)
             {
             ////cout<<i<<j<<"newStateProb"<<stateprob.first<<","<<stateprob.second<<endl;
             //if(newStateProb.find(k))
             //stateprob.second*=flowProb[j];
             // newStateProb[k]*=flowProb[j];
             //#pragma critical
             {
             //      ////////cout<<"check"<<"="<<stateprob.first<<","<<stateprob.second<<"i="<<i<<endl;
             // it_in_trans=Trans.InTrans(stateprob.first,i);
             Trans(row,k)+=newStateProb(row);
             //////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
             }
             }*/
        }
    }
    Trans.prune(0,0);
    // return Trans;
}

TransProb model3::transComputation_ignore(int ignored_flow)
{
    //omp_set_num_threads(2);
    /* if(Trans.size() ==0){
     transComputation();
     }*/
    TransProb newtrans=Trans;
    //TransProb::iterator it_in_trans;
//#pragma omp parallel for
    for( StateType k=0;k<legalState.size();++k){
        // StateType i=legalState[k];
        int j=ignored_flow;
        StateProb2 newStateProb=allnewStateProb[k][j];
        //  long double baseprob=flowProb(j);
        //  ////cout<<"newStateProb"<<endl;
        
        // ////cout<<"for"<<newStateProb.size()<<endl;
        //for( StateType k=0;k<legalState.size();++k){
        //newStateProb*=baseprob;
        // newtrans=(newtrans-newStateProb);
        newtrans.col(k)-=newStateProb;
        /*  for(size_t row = 0, nRows = newStateProb.rows(); row < nRows; ++row)
         {
         ////cout<<i<<j<<"newStateProb"<<stateprob.first<<","<<stateprob.second<<endl;
         //if(newStateProb.find(k))
         //  stateprob.second*=baseprob;
         // newStateProb[k]*=flowProb[j];
         #pragma critical
         {
         //      ////////cout<<"check"<<"="<<stateprob.first<<","<<stateprob.second<<"i="<<i<<endl;
         
         //it_in_trans=newtrans.InTrans(stateprob.first,i);
         //long double now=newtrans[it_in_trans->first];
         //////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
         if (newtrans(row,k))
         {
         //////cout<<"set"<<stateprob.first<<","<<i<<"="<<endl;
         newtrans(row,k)-=newStateProb(row);
         //it_in_trans->second -=stateprob(i);
         //  ////////cout<<"it_in_trans"<<"="<<it_in_trans->second<<endl;
         }
         }
         }*/
    }
    return newtrans;
}
void model3::init_triggerFlowP(){
    
    //   //cout<<"baseFullState"<<baseFullState<<legalState[baseFullState]<<endl;;
    triggerFlowPTable.resize((stateNum-baseFullState)*nRule);
    // StateType fullNum=(stateNum-baseFullState);
    // //cout<<"init trigger flow"<<(stateNum-baseFullState)*nRule<<endl;
    //#pragma omp parallel for
    for (StateType state=baseFullState; state<stateNum; state++) {
        StateType base=nRule*(state-baseFullState);
        StateType list=legalState[state];
        // //cout<<"init trigger flow"<<state<<endl;
        for (int i=1; i<=nRule; ++i) {
            bool exist=exist_bit(list,i);
            //long double tmp=triggerFlowP( i, list,exist)+(exist?0:MASK_EXIST);
            triggerFlowPTable[base+i-1]=triggerFlowP( i, list,exist);
        }
    }
    //   //cout<<"init trigger flow"<<endl;
    // //cout<<triggerFlowPTable<<endl;
    // triggerFlowPTable.reserve();
    //triggerFlowPTable.setZero();
}


long double model3::TTLProb_reuse(int rule, StateType stateno)
{
    StateType state=legalState[stateno];
    int ruleNum=nRule;
    StateType triggerFlowPTablebase=(stateno-baseFullState)*nRule;
    long double prob,p;
    long double total,lambda;
    long double ttlrule=TTL->get(rule);
    //int statesize=nonZeroNum(state);
    bool full=(stateno>=baseFullState);
    if (interval < ttlrule)
    {
        prob = 0;
        return prob;
    }
    prob = 1;
    // long double *lambdas=(long double *)malloc(sizeof(long double)*ruleNum);
    vector<long double> lambdas(ruleNum);
    //memset(lambdas,1,sizeof(long double)*ruleNum);
    //   ////cout<<"lambdas"<<endl;
    // BoolState exists=(ruleNum);
    // total = 1;
    if (full){
        total = 0;
        for(int i=1; i<=ruleNum; ++i)
        {
            bool existi=exist_bit(state,i);
            //exists[i-1]=existi;
            ////////cout<<"ttlprob:for"<<endl;
            lambda =triggerFlowPTable[triggerFlowPTablebase+i-1];// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            lambdas[i-1]=lambda;
            // //cout<<"lambda"<<lambda<<endl;
            /*   if (lambda==0) {
             return 0;
             }*/
            ////cout<<lambdas[i-1]<<" "<<lambda<<endl;
            ////////cout<<"ttlprob:after triggerFlowP"<<endl;
            if (i == rule){
                prob *=(1 - poissonNumber0(lambda, 0, unit)) * poissonNumber0(lambda, 0, ttlrule - unit);
                //     //cout<<"prob"<<prob<<endl;
            }
            else {
                if (existi){
                    prob *=(1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i))));
                    // //cout<<"prob"<<prob<<endl;
                }
                else{
                    prob *= poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                    //  //cout<<"prob"<<prob<<endl;
                }
            }
            /*if(prob==0)
             return 0;*/
        }
        //   //cout<<"prob"<<prob<<endl;
        
        int it;
        //  next_valid_rule(state>>=1,state);
        total=0;
        forAllRinS(it, state)
        {
            int ceil0=ceilM(min(ttlrule,TTL->get(it)), unit, delta);
            for (int k =mSize; k<=ceil0; ++k)
            {
                p = 1;
                for(int j=1; j<=ruleNum; ++j)
                {
                    lambda=lambdas[j-1];
                    if (!exist_bit(state,j))
                        p*=poissonNumber0(lambda, 0, min(k * unit, TTL->get(j)));
                    else
                    {
                        if (j == it)
                            p*=poissonNumber0(lambda, 0, (k - 1) * unit) * (1 - poissonNumber0(lambda, 0, unit));
                        else
                            p *=(1 - poissonNumber0(lambda, 0, min(k * unit, TTL->get(j))));
                    }
                    //////cout<<ceil0<<"k="<<"lambda"<<lambda<<"p="<<p<<endl;
                }
                
                total+= p;
            }
        }
        ////cout<<"lambdas3"<<endl;
    }else{
        total = 1;
        for(int i=1; i<=ruleNum; ++i)
        {
            bool existi=exist_bit(state,i);
            ////////cout<<"ttlprob:for"<<endl;
            // ////cout<<"lambdas2"<<i<<state<<existi<<endl;
            lambda =triggerFlowP(i,state,existi);// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
            //  ////cout<<"lambdas4"<<endl;
            lambdas[i-1]=lambda;
            /* if (lambda==0) {
             return 0;
             }*/
            ////////cout<<"ttlprob:after triggerFlowP"<<endl;
            if (i == rule){
                prob*=(1 - poissonNumber0(lambda, 0, unit)) * poissonNumber0(lambda, 0, ttlrule - unit);
                if (existi)
                    total*=(1 - poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i))));
                else
                    total*=poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                
            } else {
                long double tmp=poissonNumber0(lambda, 0, min(ttlrule, TTL->get(i)));
                if (existi){
                    prob*=tmp;
                    total*=tmp;
                }
                else{
                    prob*=tmp;
                    total*=tmp;
                }
            }
        }
        
    }
    //////cout<<"total="<<total<<endl;
    //////cout<<"prob="<<prob<<endl;
    ////cout<<"prob="<<prob<<endl;
    if(prob==0)
        return 0;
    if (total==0){
        for(int k=0;k<ruleNum;++k)
            //cout<<lambdas[k]<<" ";
        cerr<<"rule="<<ruleNum<<"full="<<full<<"total=="<<total<<"prob="<<prob<<endl;
        prob=0;
    }else
        prob = prob / total;
    // free(lambdas);
    return prob;
}

long double model3::ruleEVT_reuse(int rule, StateType state)
{
    StateType triggerFlowPTablebase=(state-baseFullState)*nRule;
    //////cout<<"ruleEVT"<<endl;
    int ruleNum=flowRuleTable->get_rulenum();
    long double ttlrule=TTL->get(rule);
    long double prob = 0,p;
    StateType list=legalState[state];
    /* if (list.size() > mSize)
     cerr<< "larger than cache size!";*/
    prob = 0;
    int maxv=ceilM(ttlrule, unit, delta);
    for (int k = mSize; k<=maxv; ++k)
    {
        p = 1;
        for(int i = 1; i<=ruleNum; ++i)
        {
            long double rPara =triggerFlowPTable[triggerFlowPTablebase+i-1];
            bool exist=exist_bit(list,i);
            if (!exist)
            {
                p *=poissonNumber0(rPara, 0, min(k*unit, TTL->get(i)));
            }
            else
            {
                if (rule == i){
                    p *= poissonNumber0(rPara, 0, (k -1)* unit) * poissonNumber1(rPara, 1, unit);
                }
                else{
                    p *= (1 - poissonNumber0(rPara, 0, min(TTL->get(i), k*unit)));
                }
            }
        }
        prob+=p;
    }
    //////cout<<"ruleEVT"<<endl;
    return prob;
}
void model3::flowState(int flow, StateType oldStateNum,StateProb2 & newStateProb)
{
    //    ////////cout<<"flowState3:enter"<<endl;
    StateType oldbinState=legalState[oldStateNum];
    bool needEvict=((oldStateNum>= baseFullState)&&(flow > 0));
    int ruleNum=nRule;
    LISTINT oldList,midList;
    int matched_prio=0;//matched_rule=0;
    int matched_rule=0;
    int high_prio_rule=flowRuleTable->get_high_rule(flow);
    StateProb2 midStateProb(stateNum),ttlStateProb(stateNum);
    midStateProb.setZero();
    ttlStateProb.setZero();
    newStateProb.setZero();
    size_t ttl_firststate=-1;
    matched_rule=doesMatch(flow,oldbinState);
    if (needEvict)
    {
        if(!matched_rule){
            //////cout<<"flowState3:if"<<endl;
            //LISTINT::iterator it=oldList.begin();
            int deletedrule;
            // StateType total=oldbinState;
            long double totalmid=0;
            //#pragma omp parallel for
            for (int deletedrule=1; deletedrule<=nRule; deletedrule++) {
                if(exist_bit(oldbinState, deletedrule)){
                    StateType oneBinNum=clear_bit(oldbinState,deletedrule);
                    StateType newBinNum=cacheLRU(flow, oneBinNum,(matched_rule!=deletedrule)?matched_rule:0);
                    
                    long double tmp=ruleEVT(deletedrule, oldbinState,true);
                    
                    // long double tmp=ruleEVT_reuse(deletedrule, oldStateNum);
                    ////cout<<"tmp="<<tmp<<"deltedt"<<deletedrule<<endl;
                    //#pragma omp critical
                    {
                        midStateProb.insert(bin2Num(newBinNum))=tmp;
                    }
                }
            }
            midStateProb/=midStateProb.sum();
        }else{
            StateType newBinNum=cacheLRU(flow, oldbinState ,matched_rule);
            midStateProb.setZero();
            midStateProb.insert(bin2Num(newBinNum))=1;
        }
        // ////////cout<<"flowState3:norm"<<endl;
    }
    else
    {
        StateType newBinNum;
        newBinNum=cacheLRU(flow, oldbinState,matched_rule);
        //    ////cout<<"flowState3:cache2"<<endl;
        midStateProb.setZero();
        midStateProb.insert(bin2Num(newBinNum))=1;
    }
    
    // midStateProb.prune(epsilon);
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
        if (it.index()>=baseNum[k-1])
        {
            //  ////cout<<"noless"<<endl;
            vector<long double> lambdas(nRule);
            for(int i=1; i<=nRule; ++i)
            {
                bool existi=exist_bit(newList,i);
                //triggerFlowPTable(i,it.index())
                long double lambda0 =triggerFlowP(i,newList,existi);// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
                lambdas[i-1]=lambda0;
            }
            
            
            long double ptotal=TTLStateProb(newList,lambdas);
            
            ttlStateProb.insert(it.index())=it.value();
            ttl_firststate=it.index();
            int j_it;
            forAllRinS(j_it, oldbinState) {
                //if((j_it!=newrule)&&exist_bit(oldbinState, j_it)){
                // //cout<<"ttl table"<<j_it<<"index="<<it.index()<<ttlProbTable;
                if(j_it==newrule){
                    continue;
                }
                // long double p =TTLProb_reuse(j_it, it.index());
                
                long double p =TTLProb(j_it, newList,lambdas);
                
                /*     long double p=ttlProbTable.coeffRef(it.index(),j_it-1);
                 if(p){
                 p=(p>MASK_EXIST)?0:p;
                 }else{
                 //p =TTLProb(j_it, newList);
                 p =TTLProb_reuse(j_it, it.index());
                 ttlProbTable.insert(it.index(),j_it-1)=(p?p:MASK_EXIST);
                 }*/
                if (p > epsilon)
                {
                    StateType onelist=clear_bit(newList, j_it);
                    //cout<<p<<endl;
                    p/=ptotal;
                    long double pp=p*it.value();
                    if (pp>1) {
                        pp=1;
                        caseType=CASE_SPECIAL;
                    }
                    ////cout<<p<<"newlist"<<newList<<"total="<<ptotal<<"pp="<<pp<<endl;
                    ttlStateProb.insert(bin2Num(onelist))=pp;
                    ttlStateProb.coeffRef(ttl_firststate)-=pp;
                    // ttlStateProb.coeffRef(ttl_firststate)-=ttlStateProb.sum();
                }
                
            }
            
            
        }else
        {
            ttlStateProb.coeffRef(it.index())=it.value();
        }
        newStateProb+=ttlStateProb;
    }
    newStateProb.prune(0);
    // //cout<<"end flowstate"<<endl;
    // newStateProb.prune(epsilon);
}

void model::initFlowProb(){
    
    flowProb.resize(nFlow+1);
    flowProb.setZero();
    long double prob=1;
  //  //cout<<"unit="<<unit;
    for (int i = 1; i<=nFlow; ++i)
    {
        prob *= poissonNumber0(flowPara->get(i), (int)0,unit);
        cout<<"prob"<<prob<<endl;
    }
    
    ////////cout<<"hello"<<endl;
    flowProb(0)=prob;
    //#pragma omp parallel for
    for(int i = 1; i<=nFlow; ++i)
    {
        
        flowProb(i) = poissonNumber1(flowPara->get(i), (int)1, unit) * prob/ poissonNumber0(flowPara->get(i), (int)0, unit);
    }
    flowProb/=flowProb.sum();
    //probNormalization(flowProb,1+nFlow);
}


void model3::run()
{
    stateProb.resize(stateNum,1);
    stateProb.reserve(stateNum);
    stateProb.setZero();
    stateProb.insert(initialStateNum) = 1;
    
    //deque<long double> flowProb(1+nFlow,1);
    // = (long double*)malloc((1+flowNum)*sizeof(long double));
    
    
    ////cout<<"build flowProb finished"<<endl;
    ////////cout<<"hello"<<endl;
    transComputation();
    ////////cout<<"hello end trans"<<endl;
    // long double sumtrans[Trans.size()];
    // memset(sumtrans)
    ////cout<<"build trans finished"<<endl;
    for (int i=0; i<fn; ++i) {
        stateProb=Trans*stateProb;
    }
    
    /*for(int i = 1;i<=fn;++i){
     //         ////////cout<<"matMultiply"<<endl;
     matMultiply(Trans, stateProb);
     }*/
    //  ////////cout<<"matMultiply end"<<endl;
    ////////cout<<"hello statenhum="<<stateNum<<endl;
    //return stateProb;
}
