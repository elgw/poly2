/* Lua interface, should be compiled as .so object */
/*

lpoly = require "lpoly"

L = {1, 2, 3, 4, 5, 6, 7, 8}

i = lpoly.lines_intersect(L)

lpoly.lines_intersect(1)

*/

#include "poly.h"
// #define lua51
// #define lua53

#ifdef lua51
#include "lua5.1/lua.h"
#include "lua5.1/lauxlib.h"
#define lua_version_set
#endif

#ifdef lua53
#include "lua5.3/lua.h"
#include "lua5.3/lauxlib.h"
#define lua_version_set
#endif

#ifndef lua_version_set
#error Neither lua51 or lua53 was set
#include <stophere>
#endif

static int l_lines_intersect (lua_State *L) {

    /* We expect a table with 8 elements */
    luaL_checktype(L, 1, LUA_TTABLE);
#ifdef lua51
    int npoints = lua_objlen(L, 1);
#else
    int npoints = luaL_len(L, 1);  /* get size of table */
#endif

    if (npoints != 8) { /* Wrong number of points */
        lua_pushnil(L); /* return nil... */
        lua_pushstring(L, "Wrong number of points, requires 4, i.e. 8 coordinates");
        return 2; /* number of results */
    }

    /* Get the values */
    double coords[8];
    // lua_isnumber()
    for(int kk = 1; kk<9; kk++)
    {
        // Push on the stack
        lua_rawgeti (L, 1, kk);
        double data = lua_tonumber(L, -1);
        //printf("%f\n", data);
        coords[kk-1] = data;
        // Need to pop?
    }

    int ans = lines_intersect(coords, coords+2, coords+4, coords+6);
    lua_pushnumber (L, ans);
    return 1;
}

static int l_poly_is_simple(lua_State *L) {

        /* We expect a table with 8 elements */
        luaL_checktype(L, 1, LUA_TTABLE);
#ifdef lua51
        int npoints = lua_objlen(L, 1);
#else
        int npoints = luaL_len(L, 1);  /* get size of table */
#endif

        if (npoints < 6) { /* Wrong number of points */
            lua_pushnil(L); /* return nil... */
            lua_pushstring(L, "Too few points, requires at least 3");
            return 2; /* number of results */
        }
        if (npoints % 2 != 0) {
            lua_pushnil(L); /* return nil... */
            lua_pushstring(L, "Need an even number of coordinates");
            return 2; /* number of results */
        }

        /* Get the values */
        double * coords = malloc(npoints*sizeof(double));
        // lua_isnumber()
        for(int kk = 1; kk<=npoints; kk++)
        {
            // Push on the stack
            lua_rawgeti (L, 1, kk);
            double data = lua_tonumber(L, -1);
            //printf("%f\n", data);
            coords[kk-1] = data;
            // Need to pop?
        }

        int ans = poly_is_simple(coords, npoints/2);
        free(coords);
        lua_pushnumber (L, ans);
        return 1;
    }


static const struct luaL_Reg lpoly [] = {
    {"lines_intersect", l_lines_intersect},
    {"poly_is_simple", l_poly_is_simple},
    {NULL, NULL}};

int luaopen_lpoly(lua_State *L) {
#ifdef lua51
    luaL_register(L, "lpoly", lpoly);
#else
    luaL_newlib(L, lpoly);
#endif
    return 1;
}
