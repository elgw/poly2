-- For the Löve game engine, run with
-- & love .
--
-- Note that love is built using an old version of lua, see the _VERSION string

lpoly = require "lpoly"

print(_VERSION)

w = 600 -- width of screen

points = {}
dragging = 0
updatepos = 0

printpoints = function(tab)
   -- print(type(tab))
   if (type(tab) == 'number') then
      io.write(tostring(" " .. tab))
      return
   end
   if (type(tab) == 'table') then
      io.write("[")
      for k, v in pairs(tab) do
         -- print("k=" .. k)
         printpoints(tab[k])
      end
      io.write("]\n")
   end
   end

drawhull = function(hull)
   printpoints(hull)
   love.graphics.setColor(.5,.5,.5)
   npoints = #hull
   print("Hull size: " .. tostring(npoints))
   if npoints > 2 then
      for i=1,npoints-1 do
         p1 = hull[i]
         p2 = hull[i+1]
         line = {p1[1], p1[2], p2[1], p2[2]}
         love.graphics.line(line)
      end
      p1 = hull[npoints]
      p2 = hull[1]
      line = {p1[1], p1[2], p2[1], p2[2]}
      love.graphics.line(line)
   end
end

love.draw = function()



   simple = 0
   if lpoly.poly_is_simple(points) == 1 then
      simple = 1
   else
      love.graphics.print(tostring(#points/2) .. " points", 10, 10)
      love.graphics.print("Simple: NO", 10, 40)
   end


   if simple==1 then
      properties = lpoly.poly_measure(points)
      love.graphics.print(properties, 10, 10)
      hull = lpoly.poly_hull(points)
      drawhull(hull)
      com = lpoly.poly_com(points)
      love.graphics.setColor(1,0,1)
      love.graphics.setPointSize(2)
      love.graphics.points(com[1], com[2])
   end

   if simple == 1 then
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
   if dragging==1 then
      points[updatepos*2-1] = love.mouse.getX()
      points[updatepos*2] = love.mouse.getY()
   end
end

function love.mousepressed(x, y, button, istouch)
   if button == 1 then -- Versions prior to 0.10.0 use the MouseConstant 'l'
      table.insert(points, x)
      table.insert(points, y)
   end
   if button == 2 then
      dragging = 1
      updatepos = 1;
      -- find closest point
      npoints = #points/2
      mindist = w+w
      for i=1,npoints do
         dx = (x-points[2*i-1])
         dy = (y-points[2*i])
         dist = dx*dx + dy*dy
         if dist<mindist then
            mindist = dist
            updatepos = i
         end
      end
      -- print("Updatepos=" .. tostring(updatepos))
   end
end

function love.mousereleased(x,y,button)
   if button == 2 then
      dragging = 0
   end
end

function love.keypressed(key)
   if key == "escape" then
      love.event.quit()
   end
end
