-- For the Löve game engine, run with
-- & love .
--
-- Note that love is built using an old version of lua, see the _VERSION string

lpoly = require "lpoly"

print(_VERSION)

w = 600 -- width of screen

points = {}

love.draw = function()

   if lpoly.poly_is_simple(points) == 1 then
      love.graphics.setColor(0,1,0)
   else
      love.graphics.setColor(1,0,0)
   end

   npoints = #points/2
   if npoints > 2 then
      for i=1,npoints-1 do
         line = {points[2*i-1], points[2*i], points[2*i+1], points[2*i+2]}
         love.graphics.line(line)
      end
      line = {points[npoints*2-1], points[npoints*2], points[1], points[2]}
      love.graphics.line(line)
   end
   love.graphics.setColor(1,1,1)
   love.graphics.setPointSize(2)
   for i = 1,npoints do
      love.graphics.points(points[2*i-1], points[2*i])
   end

   -- i = lpoly.lines_intersect(c)
end


love.update = function(dt)
   if love.mouse.isDown(2) then
      love.event.quit()
   end
end

function love.mousepressed(x, y, button, istouch)
   if button == 1 then -- Versions prior to 0.10.0 use the MouseConstant 'l'
      table.insert(points, x)
      table.insert(points, y)
   end
   if button == 2 then
      -- remove closest points or something ...
   end
end


function love.keypressed(key)
   if key == "escape" then
      love.event.quit()
   end
end
