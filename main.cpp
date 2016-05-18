#include <iostream>
#include "sdnsim.h"
#include<list>
#include<deque>
#include"model3.h"
#include "attacker.h"
#include <set>
#include <vector>
#include<cmath>
#include<bitset>
#include "autopara.h"
/*StateType clear_bit(StateType u,int b){
 StateType t=(1<<(b-1));
 t=~t;
 return u&t;
 };
 StateType set_bit(StateType u,int b){
 StateType t=(1<<(b-1));
 return u|t;
 };*/

#define entropy(pr) ((pr<=0)||(pr>=1))?0:(-pr*log(pr)-(1-pr)*log(1-pr))
using namespace std;
int main()
{
    StateProb2 test(4);
    //cout<<test;
    long double pr=0.5;
       vector<int> state=bin2SetAttack(11, 4, 3);
    for(int i=0;i<state.size();i++){
        cout<<state[i]<<" ";
    }
    cout<<endl;
     state=bin2SetAttack(14, 4, 3);
    for(int i=0;i<state.size();i++){
        cout<<state[i]<<" ";
    }
    cout<<endl;
    //cout<<"factory"<<factorial(13,12)<<endl;
    // //cout<<"entropy="<<(-pr*log(pr)-(1-pr)*log(1-pr))<<endl;
    //////cout << "Machine Epsilon is: " << numeric_limits<long double>::epsilon() << endl;
    //////cout<<nChoosek(10, 3)<<endl;
    //////cout<<poissonNumber(0.5,0,0.025)<<endl;
    //////cout<<ceilM(1,0.125,0.01)<<endl;
    int flowNum=3,ruleNum=3;
    int mSize=2;
    StateType initialStateNum=0;
    long double interval=2;
    long double limit=0.001;
    long double delta=0.001;
    int target=1;
    int qnum=1;
    FlowRuleTable table(flowNum,ruleNum);//flowNUm,ruleNUm
    floatCounter flowPara(flowNum);
    
    table<<2, 1, 0,
    0, 1, 0,
    0, 0, 1;
    
    
    //table.set(1,1,3);
    //table.set(2,2,2);
    //table.set(2,3,1);
    //table.set(3,3,1);
    //table.set(3,1,2);
    //table.set(3,3,3);
    //table.set(4,4,3);
    //table.set(5,1,2);
    //table.set(5,3,1);
    //table.set(5,4,3);
    //int tmp;
    //int rule;
    //  //cout<<table.row(1).maxCoeff(&tmp,&rule)<<"rule="<<rule;
    flowPara[1]=0.3;
    flowPara[2]=0.8;
    flowPara[3]=0.5;
   // flowPara[4]=0.003;
    //flowPara[4]=0.4;
    //flowPara[5]=0.1;
    //flowPara[4]=0.2;
    floatCounter TTL(ruleNum);
    TTL[1]=1;
    TTL[2]=1;
    TTL[3]=1;
    //TTL[4]=0.5;
    // for(int i=4;i<=ruleNum;++i){
    //     table.set(i,i,1);
    
    // }
    // for(int i=1;i<=ruleNum;++i){
    //     TTL[i]=1;
    // }
    // for(int i=4;i<=flowNum;++i)
    //     flowPara[i]=0.4;
    //table.set(3,1,1);*/
    /*
     table.set(4,4,1);
     table.set(5,5,1);
     table.set(6,6,1);
     table.set(7,7,1);
     table.set(8,8,1);
     table.set(9,9,1);
     table.set(10,10,1);*/
    // flowPara<<0.3,0.8,0.5;
    
    /*
     for(int i=4;i<=flowNum;++i)
     flowPara[i]=0.4;*/
    //cout<<table;
    //////cout<<"table"<<endl;
    //////cout<<"table"<<endl;
    // LISTINT state;
    //LISTINT::iterator it;
    //////cout<<"nonzero="<<countn<<endl;
    //////cout<<"num2state2"<<state.empty()<<state.size()<<endl;
    FlowRuleTable *flowRuleTable=&table;
    //  ////cout<<"high prio"<<flowRuleTable->get_high_rule(1);
    //  //cout<<"init TTl"<<endl;
    //TTL<<1,1,1;
    /*  for(int i=1;i<=ruleNum;++i){
     TTL[i]=1;
     }*/
    //   //cout<<"init TTl"<<endl;
    //  //cout<<TTL;
    ////cout<<*flowRuleTable;
    ////cout<<flowPara;
    //   //cout<<"init TTl"<<endl;
    long double unit=unitComputation(&flowPara, delta, limit, &TTL);
   
    //  //cout<<"init unit"<<endl;
    
    
    //  //cout<<"init unit"<<endl;
    //cout<<"unit="<<unit<<endl;
#if 0
    model3 model(&flowPara,flowRuleTable,&TTL,mSize,initialStateNum,interval,unit,delta);
    ////cout<<"end"<<endl;
    
    //model.run();
    ////cout<<"end"<<endl;
    StateProb2 stateprob=model.getStateProb();
    // int stateid=0;
    //char* c=(char *)malloc(10);
    // char* c=(char *)malloc(10);
    //memset(c,1,10);
    StateType l=1023;
#endif
    set<vector<int>>attackFlow;
    MatD PrXQ;
    VecD IG;
#if 1
    Attacker attacker(&flowPara,flowRuleTable,&TTL,mSize,initialStateNum,interval,unit,delta);
    ////cout<<"attack"<<endl;
    ////cout<<"run1"<<endl;
    attacker.run(qnum,target,initialStateNum,attackFlow,PrXQ,IG);
    cout<<IG<<endl;
    cout<<PrXQ<<endl;
    cout<<attacker.flowProb<<endl;
    cout<<attacker.unit<<endl;
    cout<<attacker.Trans;
    cout<<attacker.stateProb;
    StateProb2 stateprob=attacker.getStateProb();
    //cout<<attacker.get_nrule()<<endl;
    // //cout<<factorial(13)<<endl;;
    ////cout<<factorial(12)<<endl;;
  //  //cout<<(factorial(13)/factorial(12))<<endl;
 //   //cout<<"flownum="<<attacker.flowRuleTable->get_flownum()<<endl;
  //  //cout<<"flownum="<<attacker.queryNum<<endl;
    // //cout<<"prob state="<<endl;
    // //cout<<stateprob;
    // int stateid=0;
    /* for(int i=0;i<attacker.stateNum;++i){
     //cout<<i<<attacker.legalState[i]<<endl;
     }*/
#endif
    // Automatic a(&flowPara,flowRuleTable,&TTL,mSize,initialStateNum,interval,unit,delta);
    
    // int interestflow=a.generate();
    /*
     //cout<<attacker.total_time<<endl;
     vector<long double> lambdas(ruleNum);
     StateType newList=5;
     for(int i=1; i<=ruleNum; ++i)
     {
     bool existi=exist_bit(newList,i);
     long double lambda0 =attacker.triggerFlowP(i,newList,existi);// triggerFlowP(flowPara, i, flowRuleTable, state,bool_state);
     lambdas[i-1]=lambda0;
     
     }
     
     //cout<<"ceil"<<ceilM(1, unit, delta)<<endl;
     //cout<<attacker.TTLStateProb(5,lambdas);*/
    // ////cout<<"triger flowp"<<triggerFlowP(flowPara, 2, flowRuleTable, oldList,bool_oldList);
    //////cout<<"maxv=ceilM(TTL[rule], unit, delta)"<<ceilM(TTL[1], unit, delta)<<endl;
    return 0;
}

