//
//  lua_pmain.h
//  Mechatronix
//
//  Created by Paolo Bosetti on 9/24/13.
//  Copyright (c) 2013 UniTN. All rights reserved.
//

#ifndef __Mechatronix__lua_pmain__
#define __Mechatronix__lua_pmain__
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define lua_c

#ifdef ENABLE_READLINE
#define LUA_USE_READLINE
#endif


#include <lua.hpp>

int  pmain (lua_State *L);
void finalreport (lua_State *L, int status);
void l_message (const char *pname, const char *msg);

#endif /* defined(__Mechatronix__lua_pmain__) */
