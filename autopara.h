//
//  autopara.h
//  sdnmodel
//
//  Created by ziqiao zhou on 5/14/16.
//  Copyright (c) 2016 ziqiao zhou. All rights reserved.
//

#ifndef __sdnmodel__autopara__
#define __sdnmodel__autopara__

#include <stdio.h>
#include "attacker.h"
class Automatic:public Attacker{
    void automatic();
    void paraGenerate(int flowNum,int ruleNum,double alpha,double TTLMax);
};

typedef int[2] RulePrio;
#endif /* defined(__sdnmodel__autopara__) */
