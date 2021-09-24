-- For the Löve game engine, run with
-- & love .
--
-- Note that love is built using an old version of lua, see the _VERSION string

lpoly = require "lpoly"

print(_VERSION)

it = 0 -- integrated time since last update
w = 600 -- width of screen
i = 0 -- if they intersect
l1 = {50, 0, 0, 100}
l2 = {0, 100, 200, 200}

join_lines = function(a, b)
   c = {l1[1], l1[2], l1[3], l1[4], l2[1], l2[2], l2[3], l2[4]}
   return c
end

love.draw = function()
   c = join_lines(l1, l2)
   love.graphics.line(l1)
   love.graphics.line(l2)
   i = lpoly.lines_intersect(c)
end


love.update = function(dt)
   -- print(dt)
   it = it + dt
   if i == 1 then
      it = 0
      l1 = {math.random()*w, math.random()*w, math.random()*w, math.random()*w}
      l2 = {math.random()*w, math.random()*w, math.random()*w, math.random()*w}
      i = 0
   end
   l1[1] = love.mouse.getX()
   l1[2] = love.mouse.getY()
   if love.mouse.isDown(1) then
      -- print('love.mouse.isDown(1)')
      love.event.quit()
   end
end
