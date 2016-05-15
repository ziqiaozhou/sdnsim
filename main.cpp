#include <iostream>
#include "sdnsim.h"
#include<list>
#include<deque>
#include"model3.h"
#include "attacker.h"
#include <set>
#include <vector>
#include<cmath>
StateType clear_bit(StateType u,int b){
    StateType t=(1<<(b-1));
    t=~t;
    return u&t;
};
StateType set_bit(StateType u,int b){
    StateType t=(1<<(b-1));
    return u|t;
};

#define entropy(pr) ((pr<=0)||(pr>=1))?0:(-pr*log(pr)-(1-pr)*log(1-pr))
using namespace std;
int main()
{
    StateProb2 test(4);
    cout<<test;
    double pr=0.5;
    cout<<"factory"<<factorial(13,12)<<endl;
      // cout<<"entropy="<<(-pr*log(pr)-(1-pr)*log(1-pr))<<endl;
    ////cout << "Machine Epsilon is: " << numeric_limits<double>::epsilon() << endl;
    ////cout<<nChoosek(10, 3)<<endl;
    ////cout<<poissonNumber(0.5,0,0.025)<<endl;
    ////cout<<ceilM(1,0.125,0.01)<<endl;
    int flowNum=3,ruleNum=3;
    int mSize=2;
    StateType initialStateNum=0;
    double interval=2;
   	double limit=0.001;
    double delta=0.001;
    int target=1;
    int qnum=1;
    FlowRuleTable table(flowNum,ruleNum);//flowNUm,ruleNUm
    floatCounter flowPara(flowNum);
  //  cout<<table;
 /*   table<< 2,0,0,0,0,0,0,
            2,3,0,4,7,5,0,
            0,3,1,0,7,0,0,
            2,3,0,0,0,0,0,
            0,0,1,0,0,5,0,
            0,0,1,0,0,0,6,
            2,0,1,0,0,5,0,
            0,0,0,0,0,0,6,
            0,0,1,4,0,0,6,
            0,0,0,4,0,0,0;
    flowPara[1]=0.8979;
    flowPara[2]=0.1363;
    flowPara[3]=0.5778;
    flowPara[4]=0.3874;
    flowPara[5]=0.5930;
    flowPara[6]=0.6874;
    flowPara[7]=0.3825;
    flowPara[8]=0.0628;
    flowPara[9]=0.3527;
    flowPara[10]=0.8408;
    TTL[1]=0.4;
    TTL[2]=1;
    TTL[3]=0.2;
    TTL[4]=0.6;
    TTL[5]=0.2;
    TTL[6]=0.4;
    TTL[7]=0.6;
*/
    table.set(1,1,2);
    ////cout<<"table"<<endl;
   table.set(1,2,1);
    //table.set(4,1,1);
    //table.set(4,3,2);
    ////cout<<"table"<<endl;
    table.set(2,2,1);
    ////cout<<"table"<<endl;
    table.set(3,3,1);
    flowPara[1]=0.3;
    flowPara[2]=0.8;
    flowPara[3]=0.5;
    //flowPara[4]=0.2;
    floatCounter TTL(ruleNum);

    for(int i=4;i<=ruleNum;i++){
        table.set(i,i,1);
        
    }
    for(int i=1;i<=ruleNum;i++){
        TTL[i]=1;
    }
    for(int i=4;i<=flowNum;i++)
        flowPara[i]=0.4;
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
    for(int i=4;i<=flowNum;i++)
        flowPara[i]=0.4;*/
    cout<<table;
    ////cout<<"table"<<endl;
    ////cout<<"table"<<endl;
   // LISTINT state;
    //LISTINT::iterator it;
    ////cout<<"nonzero="<<countn<<endl;
    ////cout<<"num2state2"<<state.empty()<<state.size()<<endl;
        FlowRuleTable *flowRuleTable=&table;
  //  //cout<<"high prio"<<flowRuleTable->get_high_rule(1);
  //  cout<<"init TTl"<<endl;
       //TTL<<1,1,1;
 /*  for(int i=1;i<=ruleNum;i++){
        TTL[i]=1;
    }*/
  //   cout<<"init TTl"<<endl;
     //  cout<<TTL;
    //cout<<*flowRuleTable;
    //cout<<flowPara;
  //   cout<<"init TTl"<<endl;
    double unit=unitComputation(&flowPara, delta, limit, &TTL);
   //  cout<<"init unit"<<endl;
    cout<<"unit="<<unit<<endl;
#if 0
    model3 model(&flowPara,flowRuleTable,&TTL,mSize,initialStateNum,interval,unit,delta);
    //cout<<"end"<<endl;

   model.run();
    //cout<<"end"<<endl;
    StateProb2 stateprob=model.getStateProb();
   // int stateid=0;
    //char* c=(char *)malloc(10);
    // char* c=(char *)malloc(10);
    //memset(c,1,10);
    StateType l=1023;
#endif
    set<int>attackFlow;
    vector<vector<double>>*PrXQ=NULL;
    vector<double>*IG=NULL;
    Attacker attacker(&flowPara,flowRuleTable,&TTL,mSize,initialStateNum,interval,unit,delta);
    cout<<"attack"<<endl;
    cout<<"run1"<<endl;
    attacker.run(qnum,target,initialStateNum,attackFlow,PrXQ,IG);
    for(std::set<int>::iterator it=attackFlow.begin(); it!=attackFlow.end(); ++it){
        cout<<"choose"<<*it<<endl;
    }
    cout<<"choose end"<<endl;
    
  /*  for(int i=0;i<PrXQ->size();i++){
        for(int j=0;j<(*PrXQ)[i].size();j++){
            cout<<i<<","<<j<<"="<<(*PrXQ)[i][j]<<endl;
        }
    }*/
    StateProb2 stateprob=attacker.getStateProb();
    cout<<attacker.get_nrule()<<endl;
   // cout<<factorial(13)<<endl;;
    //cout<<factorial(12)<<endl;;
    cout<<(factorial(13)/factorial(12))<<endl;
    cout<<"flownum="<<attacker.flowRuleTable->get_flownum()<<endl;
    cout<<"flownum="<<attacker.queryNum<<endl;
    cout<<"prob state="<<endl;
    cout<<stateprob;
    // int stateid=0;
       for(int i=0;i<IG->size();i++){
        set<int> set=num2SetAttack(i,flowNum,1);
        for(auto & one: set){
        cout<<one<<endl;
        }
        cout<<"gain="<<(*IG)[i]<<endl;;
    }
    for(int i=0;i<PrXQ->size();i++){
        for(int j=0;j<(*PrXQ)[i].size();j++)
            cout<<(*PrXQ)[i][j]<<" ";
        cout<<endl;
    }
  /* for(int i=0;i<attacker.stateNum;i++){
        cout<<i<<attacker.legalState[i]<<endl;
    }*/
    free(PrXQ);
    free(IG);
 
    // //cout<<"triger flowp"<<triggerFlowP(flowPara, 2, flowRuleTable, oldList,bool_oldList);
    ////cout<<"maxv=ceilM(TTL[rule], unit, delta)"<<ceilM(TTL[1], unit, delta)<<endl;
    return 0;
}
