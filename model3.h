#ifndef MODEL3_H_INCLUDED
#define MODEL3_H_INCLUDED
#include<deque>
#include"sdnsim.h"
#include<set>
#define exist_bit( u, b) u&(1<<(b-1))
#define exist_bit_c( u, b) u[b]
class model3:public model{
    public:
    //Eigen::SparseVector<double> allnewStateProb;
    std::deque<std::deque<StateProb2>>allnewStateProb;
    StateType ttlPNum;
     TransProb ttlProbTable;
    StateType stateNumCompute(int mSize,int ruleNum)
    {
        StateType stateNum = 1;
        ttlPNum=1;
        std::cout<<"resize"<<std::endl;
        int times=(mSize<ruleNum)?mSize:ruleNum;
        for (int i= 1; i<=times; i++){
            StateType tmp=nChoosek(ruleNum, i);
            stateNum+=tmp;
            ttlPNum+=tmp*i;
        }
        std::cout<<"resize"<<std::endl;
        return stateNum;
    }
    model3(floatCounter * flowPara, FlowRuleTable *flowRuleTable, floatCounter * TTL, int mSize, int initialStateNum, double interval, double unit, double delta):model(flowPara,flowRuleTable,TTL,mSize,initialStateNum,interval,unit,delta){
        //StateType i,j=1;
        stateNum=0;
        std::cout<<"init model3"<<std::endl;

        //#pragma omp parallel for
        StateType all=1<<nRule;
        stateNum=stateNumCompute(mSize,nRule);
        legalState.resize(stateNum);
        for( StateType i=0; i<all; i++)
        {
            if(!isLegalState(i)){
                continue;
            }
            legalState[bin2Num(i)]=i;
            //stateNum++;
        }
        initFlowProb();
        ttlProbTable.resize(stateNum,nRule);
        ttlProbTable.reserve(ttlPNum);
        
    };
    StateType bin2Num(StateType state);
    void flowState(int flow, StateType oldStateNum,StateProb2& newStateProb);
    TransProb transComputation();
    TransProb transComputation_ignore(int);
    TransProb transComputation(int ignored_flow);
    void run();
    StateProb2& getStateProb(){
        return stateProb;
    }
    ~model3(){
           }
    
};
#endif // MODEL3_H_INCLUDED
