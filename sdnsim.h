#ifndef _SDNSIM_H_
#define _SDNSIM_H_
#include<list>
#include<deque>
#include <iostream>
#include<unordered_map>
#include<array>
#include<stdlib.h>
#include <cstring>
#include<vector>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
typedef int stateNumType;
typedef unsigned long long StateType;
typedef double probType;
typedef int RuleType;
typedef std::list<RuleType> LISTINT;
typedef std::list<probType> LISTFLOAT;

class Counter{
    //one bit per rule;
    char * counter;
    int ruleNum;
    int realsize;
    public:
    Counter(int ruleNum0){
        ruleNum=ruleNum0;
         realsize=(ruleNum-1)/(sizeof(char))+1;
        counter=(char* )malloc(realsize);
        memset(counter,0,realsize);
    }
    void reset(){
        memset(counter,0,(ruleNum-1)/8+1);
    }
    bool get(int index){
         int group=index/8;
         int offset=index%8;
        return ((int)(counter[group])>>offset)&1;
    }
    int size(){
        return ruleNum;
    }
    bool operator[](int idx) {
        return get(idx);
    }

    void set(int index,bool value){
        int group=index/8;
        int offset=index%8;
        if(value)
         counter[group]=(int)(counter[group])|(1<<offset);
        else
            counter[group]=(int)(counter[group])&~(1<<offset);
    }

    ~Counter(){
        free(counter);
    }
};
typedef  Counter BinState;
class BoolCounter{
    //one bit per rule;
    bool * counter;
    int ruleNum;
    int realsize;
    public:
    BoolCounter(int ruleNum0){
        ruleNum=ruleNum0;
         realsize=ruleNum*sizeof(bool);
        counter=(bool* )malloc(realsize);
        memset(counter,0,realsize);
    }
    void reset(){
        memset(counter,0,realsize);
    }
    bool get(int rule){
        return counter[rule-1];
    }
    int size(){
        return ruleNum;
    }
    bool operator[](int rule) {
        return get(rule-1);
    }

    void set(int rule,bool value){
        counter[rule-1]=value;
    }

    ~BoolCounter(){
        free(counter);
    }
};

template <class T>
//:public Eigen::Matrix<T,1,Eigen::Dynamic>
class TCounter{
    //one bit per rule;
    T * counter;
    int ruleNum;
    int realsize;
    public:
    TCounter(){
        
    }
    
    TCounter(int ruleNum0){
      //  this->resize(1,ruleNum0);
        ruleNum=ruleNum0;
         realsize=ruleNum*sizeof(T);
        //this->setZero();
        std::cout<<"init"<<std::endl;
        counter=(T* )malloc(realsize);
      //   std:://cout<<"alloc"<<counter<<std::endl;
        memset(counter,0,realsize);
    }
    /*TCounter(TCounter& t){
        counter=(T* )malloc(realsize);
        memset(counter,0,realsize);
        ruleNum=t.ruleNum;
        realsize=t.realsize;
    };*/
   /* TCounter& operator=(TCounter t){
        counter=t.counter;
        ruleNum=t.ruleNum;
        realsize=t.realsize;
        return t;
    };*/
    inline void reset()   {
         memset(counter,0,realsize);
    }
   inline  T & get(int rule)  {
        return counter[rule-1];
    }
    int size(){
        return ruleNum;
    }
    T & operator[](int rule) {
        return get(rule);
    }

    void set(int rule,T value){
        counter[rule-1]=value;
    }

    ~TCounter(){
      //  std:://cout<<"free"<<counter<<std::endl;
      //  free(counter);
    }
};

typedef  TCounter<bool> BoolState;
typedef  TCounter<stateNumType> IntState;
typedef TCounter<double> floatCounter;

class arrayHash {
public:
  std::size_t operator()(std::array<StateType, 2> const &arr) const {
    std::size_t sum(0);
    for(auto &i : arr) sum += std::hash<StateType>()(2);
    return sum;
  }
};
typedef Eigen::SparseMatrix<double> TransProb;
//class TransProb : public std::unordered_map<std::array<StateType,2>,double,arrayHash>{
/*class TransProb: public std::vector<std::vector<double>>{
//class TransProb : public std::deque<std::deque<double,arrayHash>>{
    public:
    TransProb (){
    }
    inline void set(StateType prev,StateType now,double prob){
       // std::array<StateType,2> key={prev,now};
        (*this)[prev][now]=prob;
    }
    inline double& get(StateType prev,StateType now){
      //  std::array<StateType,2> key={prev,now};
        return (*this)[prev][now];
    }
    int InTrans(StateType prev,StateType now){
        std::array<StateType,2> index={prev,now};
       ////cout<<"intrans="<<(this->find(index)!=this->end())<<std::endl;
        return (*this)[prev][now]>0;
    }
};*/
//typedef Eigen::MatrixXd::Matrix<double, Dynamic, Dynamic> MatrixXd;
//class TransProb : public Eigen::MatrixXd::Matrix<double, Dynamic, Dynamic>{
//}



class intfloat{
    public:
int state;
double prob;
intfloat(int s,double p){
    state=s;
    prob=p;
}

};
typedef Eigen::SparseVector<double> StateProb2;
//clas StateProb2s//: public  std::unordered_map<StateType,double>{
/*class StateProb2:public std::vector<double>{
    public:
    //std::unordered_map<StateType,StateType>
    void resize(int nstate){
        this->resize(nstate);
    };
    void inline push_back(StateType s,double p){
        this->at(s)=p;
    }
};*/


//::vector<vector<unsigned int>>
class FlowRuleTable: public Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>{
    //one bit per (rule,flow) pair;
   // unsigned int ** counter;
    unsigned int ruleNum;
    unsigned int flowNum;
  //  int realsize;
    public:
    FlowRuleTable(int flowNum0,int ruleNum0){
        ruleNum=ruleNum0;
        flowNum=flowNum0;
        this->resize(flowNum,ruleNum);
        this->setZero();
   //     realsize=(ruleNum)*(flowNum);

       // counter=(unsigned int **)malloc(ruleNum*sizeof(void*));
      /*  for(int i=0;i<ruleNum;++i){
                counter[i]=(unsigned int *)malloc(flowNum*sizeof(unsigned int));
                memset(counter[i],0,flowNum*sizeof(unsigned int));
        }*/
    }
     ~FlowRuleTable(){
        for(int i=0;i<ruleNum;++i){
                //std:://cout<<"free"<<i<<"ADDRESS="<<counter[i]<<std::endl;
               //free(counter[i]);
        }
        //std:://cout<<"free"<<std::endl;
        //free(counter);
    }
    inline int& get(int flow,int rule)  {
       return (*this)(flow-1,rule-1);
    }
    inline int get_rulenum()  {
       return ruleNum;
    }
    inline int get_flownum()  {
       return flowNum;
    }
    inline int get_high_rule(int flow)  {
        int rule=0;
        int prio=0;
        for(int i=1;i<=ruleNum;++i){
            if(get(flow,i)>prio){
                rule=i;
                prio=get(flow,i);
            }
        }
        return rule;
    }
    inline void set(int flow,int rule,int value)  {
        (*this)(flow-1,rule-1)=value;
    }

};

class SeqFlow{
    double seq;
    int flow;
public:
    SeqFlow(double s,int f){
        seq=s;
        flow=f;
    }
     double get_seq() const{
        return seq;
    }
     int get_flow(){
        return flow;
    }

};

struct SeqFlowSorter{
 bool operator ()(const SeqFlow &a, const SeqFlow &b) {
        return a.get_seq()<b.get_seq();
    }
};

template <class T>
inline void maxmValueFind(std::list<T> all,T & maxm, typename std::list<T>::iterator& position)  {
    maxm=-1;
    if (!all.empty())
    {
        position = all.begin();
        maxm=*position;
        for(typename std::list<T>::iterator it=all.begin();it!=all.end();++it){
                T one=*it;
                if(one>maxm){
                    maxm=one;
                    position=it;
                }
        }
    }
};
template <class T>
inline T distance(std::list<T> vec1, std::list<T>vec2)  {
    int l_1 = vec1.size();
    int l_2 = vec2.size();
    T dist=0;
    if (l_1 != l_2)
        dist = -1;
    else{
        typename std::list<T>::iterator it1;
        typename std::list<T>::iterator it2;
        for(int i=0;i<l_1;++i){
            ++it1;
            ++it2;
            dist = dist + abs((*it1) - (*it2));
        }
    }
    return dist;
};
template <class T>
 inline void minmValueFind(std::list<T> &all,T &minm, typename std::list<T>::iterator &position)  {
    minm=-1;
    position=all.end();
    if (!all.empty())
    {
        position = all.begin();
        minm=*position;
        for( typename std::list<T>::iterator it=all.begin();it!=all.end();++it){
                T one=*it;
                if(one<minm){
                    minm=one;
                    position=it;
                }
        }
    }
};
template <class T>
inline bool ifContain(T element, std::list<T> all)  {
// determine whether the list contains a specific rule
    return std::find(all.begin(), all.end(), element) != all.end();
};

//bool ifInTrans(int prev,int now, TransProb& mat,TransProb::iterator &it);
//void combineSeq(std::list<SeqFlow>& seqflow,std::list<double> seq2,int flow);
inline int ceilM(double TTL,double unit,double delta)  
{
    return ceil(TTL/unit);
}
class model{
    public:
    StateProb2 stateProb;
    floatCounter * flowPara;
    FlowRuleTable* flowRuleTable;
    floatCounter * TTL;
    Eigen::VectorXd flowProb;
    int mSize;
    int nRule;
    int nFlow;
    int initialStateNum;
    double interval;
    double unit;
    double delta;
    int fn;
    TransProb Trans;
    std::deque<StateType>legalState;
    StateType stateNum;
    
    int get_nrule(){
        return nRule;
    }
    inline StateType cacheLRU(int flow, StateType  oldList,int match_pos)  
    {
        StateType newList=0;
        if(flow==0){
            newList=oldList;
            return newList;
        }
        if (match_pos)
        {
            //  int newRule = match_pos;
            newList=oldList;//clear_bit(oldList,match_pos);
        }
        else
        {
            int newRule=flowRuleTable->get_high_rule(flow);
            newList=set_bit(oldList,newRule);
        }
        //    if(newList.size()>mSize)
        //      cerr<<"over cache size in cacheLRU()"<<endl;
        return newList;
    }

    
    inline bool  isLegalState(StateType stateNum)  {
        if((stateNum>>mSize)==0)
            return true;
        int count=0;
        while(stateNum>0){
            count+=(stateNum&1);
            stateNum>>=1;
            if(count>mSize)
                return false;
        }
        return true;
    }
    inline StateType  clear_bit(StateType u,int b)  {
        return u&(~(1<<(b-1)));
    };
    inline StateType set_bit(StateType u,int b)   {
        //StateType t=(1<<(b-1));
        return u|(1<<(b-1));
    };
    inline bool exist_bit(StateType u,int b)  {
        return (u&(1<<(b-1)));
    }
    double ruleEVT(int rule, StateType list,bool full);
    inline int  nonZeroNum(StateType stateNum)  
    {
        int count=0;
        while(stateNum){
            count+=(1&stateNum);
            stateNum=stateNum>>1;
        }
        return count;
    };
    //void matMultiply(TransProb& mat, StateProb2& oldprob);
    inline bool nolessbit(unsigned long state,int n)  {
        int count=0;
        if ((state>>(n-1))==0) {
            return false;
        }
        while (state) {
            count+=(state&1);
            state>>=1;
            if(count>=n)
                return true;
        }
        return false;
    }
    void initFlowProb();
    int doesMatch(int flow,StateType stateNum);
       // virtual int state2Num(LISTINT &state, int ruleNum,int  mSize);
     StateType stateNumCompute(int ruleNum,int mSize);
    double triggerFlowP(int rule,StateType state,bool exit);
    double TTLProb(int rule, StateType state);
   // virtual void num2State(StateType num, int ruleNum, int mSize,LISTINT &state,BoolState& cmp);
       model(floatCounter * flowPara0, FlowRuleTable *flowRuleTable0,
                                   floatCounter * TTL0, int mSize0, int initialStateNum0, double interval0, double unit0, double delta0){
         std::cout<<"init model3"<<std::endl;
           flowPara=flowPara0;
         flowRuleTable=flowRuleTable0;
         TTL=TTL0;
         mSize=mSize0;
         initialStateNum=initialStateNum0;
         interval=interval0;
         unit=unit0;
         delta=delta0;
         nRule=flowRuleTable->get_rulenum();
         nFlow=flowRuleTable->get_flownum();
        fn = ceilM(interval, unit, delta);
                             

    };
   
  /*  template<typename FUNCTION>
    void forAllRinS(StateType numState, FUNCTION func){
            int ruleNo=0;
            while(numState){
                ruleNo++;
                if(numState&1){
                    if(!func(ruleNo))
                        return;
                }
                numState=numState>>1;
            }

    };*/
};
//int next_valid_rule(unsigned long &oldState,unsigned long &numState);
#define forAllRinS(ruleNo,numState0) StateType numState=numState0;for(ruleNo=next_valid_rule(numState,numState);numState>0;ruleNo+=next_valid_rule(numState>>=1,numState))

inline StateType  nChoosek(int n, int k)  {
    if (k == 0)
        return 1;
    return (StateType)(n * nChoosek(n - 1, k - 1)) / k;

}
 static inline StateType factorial(int x)
{
    return (x <=1?1:x*factorial(x-1));
};
 static inline StateType factorial(int x,int y)
{
    if(y>=x)
        return 1;
    return (x<=y?1:x*factorial(x-1,y));
};
inline int next_valid_rule(StateType oldState,StateType &numState)  {
    int ruleNo=1;
   // StateType oldState=oldState0;
    if(oldState==0){
        numState=oldState;
        return -1;
    }
    while(!(oldState&1)){
        if (oldState==0) {
            return 0;
        }
        oldState>>=1;
        ruleNo++;
    }
    numState=oldState;
    // std:://cout<<"num"<<numState<<std::endl;
    return ruleNo;
};

Counter list2Counter(LISTINT all,int ruleNum);

int nonZeroNum(LISTINT state);
//int  num2bin(StateType stateNumber, int ruleNumber,BoolState * state,int& nonZeroNum);
//double poissonNumber(double lambda, int number,double interval);
//void cacheLRU(int flow, LISTINT& oldList,LISTINT& newList, int mSize, FlowRuleTable* flowRuleTable,LISTINT::iterator& match_pos);
//double triggerFlowP(floatCounter& flowPara,int rule,FlowRuleTable* flowRuleTable, LISTINT& state,BoolState & bool_state);
//double ruleEVT(int rule, LISTINT list,BoolState &bool_list,FlowRuleTable* flowRuleTable,
  //            int mSize, floatCounter& flowPara, floatCounter&TTL, double unit,double delta);
//void probNormalization(StateProb2& midStateProb);
//void matMultiply(TransProb&  mat, StateProb2& oldprob);
//double TTLProb(int rule, LISTINT state,BoolState & bool_state,FlowRuleTable* flowRuleTable,
            //  floatCounter & TTL,int mSize,floatCounter& flowPara,double unit,double interval,double delta);
double unitComputation(floatCounter *flowPara,double delta,double limit, floatCounter *TTL);
inline double poissonNumber(double lambda, int number,double interval)  
{
    return (exp(-lambda*interval)*pow(lambda*interval,number)/factorial(number));
};
inline double poissonNumber0(double lambda, int number,double interval)  
{
    return (exp(-lambda*interval));
};
inline double poissonNumber1(double lambda, int number,double interval)  
{
    return (exp(-lambda*interval)*lambda*interval);
};
#endif // SDNSIM_H
