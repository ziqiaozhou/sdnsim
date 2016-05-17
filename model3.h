#ifndef MODEL3_H_INCLUDED
#define MODEL3_H_INCLUDED
#include<deque>
#include"sdnsim.h"
#include<set>
#include <chrono>
#define exist_bit( u, b) (u&(1<<(b-1)))
#define exist_bit_c( u, b) u[b]
#define MASK_EXIST 2
class model3:public model{
    public:
    //Eigen::SparseVector<double> allnewStateProb;
    std::deque<std::deque<StateProb2>>allnewStateProb;
    StateType ttlPNum;
    StateType baseFullState;
    StateType maxTrans;
    StateType baseNum[2];
     TransProb ttlProbTable;
    StateType fullNum;
    long total_time;
    std::vector<double> triggerFlowPTable;
    inline StateType stateNumCompute(int mSize,int ruleNum)
    {
        StateType stateNum = 1;
        ttlPNum=1;
        maxTrans=0;
        std::cout<<"resize"<<std::endl;
        int times=(mSize<ruleNum)?mSize:ruleNum;
        baseNum[0]=nChoosek(ruleNum,0);
        baseNum[1]=baseNum[0]+nChoosek(ruleNum,1);
        for (int i= 1; i<=times; ++i){
            StateType tmp=nChoosek(ruleNum, i);
            stateNum+=tmp;
            ttlPNum+=tmp*i;
            maxTrans+=tmp*ruleNum*nChoosek(i, 2);
        }
        baseFullState=stateNum-nChoosek(nRule,mSize);
        //cout<<"base full"<<baseFullState;
        fullNum=(stateNum-baseFullState);
        return stateNum;
    }
    double TTLStateProb( StateType state,std::vector<double> lambdas);
    void init_triggerFlowP();
    double TTLProb_reuse(int rule, StateType state);
    double epsilon;
    model3(floatCounter * flowPara, FlowRuleTable *flowRuleTable, floatCounter * TTL, int mSize, int initialStateNum, double interval, double unit, double delta):model(flowPara,flowRuleTable,TTL,mSize,initialStateNum,interval,unit,delta){
        //StateType i,j=1;
        
    };
    StateType bin2Num(StateType state);
    void flowState(int flow, StateType oldStateNum,StateProb2& newStateProb);
    void transComputation();
    TransProb transComputation_ignore(int);
    TransProb transComputation(int ignored_flow);
    double ruleEVT_reuse(int rule, StateType stateNum);
    void run();
    void init(){
        unit=unitComputation(flowPara, delta, limit, TTL);
        stateNum=0;
        total_time=0;
        std::cout<<"init model3"<<std::endl;
        fn = ceilM(interval, unit, delta);
        //#pragma omp parallel for
        StateType all=1<<nRule;
        stateNum=stateNumCompute(mSize,nRule);
        legalState.resize(stateNum);
        for( StateType i=0; i<all; ++i)
        {
            if(!isLegalState(i)){
                continue;
            }
            legalState[bin2Num(i)]=i;
            //stateNum++;
        }
        //epsilon=1/(2*stateNum);
        epsilon=0;
        initFlowProb();
        // maxTrans=nChoosek(nRule, nRule/2);
        // maxTrans
      //  ttlProbTable.resize(stateNum,nRule);
        //ttlProbTable.reserve(ttlPNum);
        //ttlProbTable.setZero();
       // init_triggerFlowP();
        std::cout<<"maxTrans"<<maxTrans<<std::endl;
    };
    StateProb2& getStateProb(){
        return stateProb;
    }
    ~model3(){
           }
    
};
#endif // MODEL3_H_INCLUDED
