#include "cddefines.h"
#include "cddrive.h"
extern "C" {
#include "lua.h"
#include "lauxlib.h"
}

namespace
{
	int lua_cdInit (lua_State *)
	{
		cdInit();
		return 0;
	}
	int lua_cdRead (lua_State *L)
	{
		size_t len;
		const char *s = luaL_checklstring(L,1,&len);
		cdRead(s);
		return 0;
	}
	int lua_cdDrive (lua_State *L)
	{
		int lgOK = cdDrive();
		lua_pushboolean(L, 0==lgOK);
		return 1;
	}
	int lua_cdOutput (lua_State *L)
	{
		size_t lenfile, lenmode;
		const char *file = luaL_checklstring(L,1,&lenfile);
		const char *mode = luaL_checklstring(L,2,&lenmode);
		cdOutput(file,mode);
		return 0;
	}
	const struct luaL_Reg lua_cloudy[] =
	{
		{"init",lua_cdInit},
		{"read",lua_cdRead},
		{"drive",lua_cdDrive},
		{"output",lua_cdOutput},
		{NULL, NULL} // Sentinel
	};
}


extern "C" int luaopen_cloudy (lua_State *L)
{
	luaL_register(L, "cloudy", lua_cloudy);
	return 1;
}

