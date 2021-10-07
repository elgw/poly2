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

static int l_poly_chull(lua_State *L) {

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

        int nhull = 0;
        double * hull = poly_hull(coords, npoints/2, &nhull);
        free(coords);

        lua_newtable(L);

        for (int i = 0; i < nhull; i++) {
            lua_newtable(L);
            lua_pushnumber(L, hull[2*i]);
            lua_rawseti(L, -2, 1);
            lua_pushnumber(L, hull[2*i+1]);
            lua_rawseti(L, -2, 2);

            lua_rawseti(L, -2, i+1);
        }

        if(hull != NULL)
        {
            free(hull);
        }
        return 1;
    }

static int l_poly_com(lua_State *L) {

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


        double * com = poly_com(coords, npoints/2);
        free(coords);
        if(com == NULL)
            return 0;

            lua_newtable(L);
            lua_pushnumber(L, com[0]);
            lua_rawseti(L, -2, 1);
            lua_pushnumber(L, com[1]);
            lua_rawseti(L, -2, 2);

            free(com);
        return 1;
    }


static int l_poly_measure(lua_State *L) {

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

        poly_props * props = poly_measure(coords, npoints/2);
        free(coords);
        char * bp;
        size_t size;
        FILE * stream = open_memstream(&bp, &size);
        poly_props_print(stream, props);
        fflush(stream);
        fclose(stream);
        lua_pushstring (L, bp);
        poly_props_free(&props);
        return 1;
    }


/* Return a table with all measurements */
static int l_poly_measure2(lua_State *L) {

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

        poly_props * props = poly_measure(coords, npoints/2);
        free(coords);
        // TODO: construct table from all measurements and return

        poly_props_free(&props);
        return 1;
    }


static const struct luaL_Reg lpoly [] = {
    {"lines_intersect", l_lines_intersect},
    {"poly_is_simple", l_poly_is_simple},
    {"poly_measure", l_poly_measure},
    {"poly_measure2", l_poly_measure2},
    {"poly_hull", l_poly_chull},
    {"poly_com", l_poly_com},
    {NULL, NULL}};

int luaopen_lpoly(lua_State *L) {
#ifdef lua51
    luaL_register(L, "lpoly", lpoly);
#else
    luaL_newlib(L, lpoly);
#endif
    return 1;
}
