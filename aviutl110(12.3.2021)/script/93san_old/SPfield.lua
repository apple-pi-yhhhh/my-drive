
--[[

特定の位置から球状に押し出すフィールドを作ります。
基本的に相対座標を操作します。
スクリプトで使用する場合はrequireで呼び出して下さい。

一度他の変数に代入するか、
local Field = require("SPfield")
Field(pos,layer,radius,strength,ret)

もしくは直接引数を渡して関数を実行します。
require("SPfield")(pos,layer,radius,strength,ret)


SphericalFieldスクリプトとは少し挙動が違います。
------引数の説明-------

SPfield([pos,layer,radius,strength,ret])

pos
  押し出したいオブジェクトの座標を{x,y,z} で指定します。
  nilの場合は現在のオブジェクトの相対座標が入ります。

layer (省略時={0,0,0})
  押し出しフィールドの位置をレイヤーで指定します。
  テーブル{x,y,z}で直接座標を指定することもできます。
  nilだと{0,0,0} になります。

radius (省略時=200)
  押し出しフィールドのサイズを指定します。
  nilだとlayerで指定したオブジェクトの拡大率を使用します。
  マイナスにすると吸着します。

strength(省略時=100)
  影響度を％で指定します。

container(省略時=0)
  1の場合押し込みフィールドになります。

ret (省略時=false)
  0以外の何らかの値を入れると戻り値のみになります。
  obj.drawで大量に描画する場合に使用します。
  x,y,z,単位化された距離, の四つを返します。
  ※引数省略時であっても戻り値を返します。

]]




local Field=function(pos,layer,radius,strength,container,ret)

  local r,f,L = radius ,strength or 100,layer
  f = f*.01
  --f=math.max(-1,math.min(2, f))
  local gv=obj.getvalue
  local pos=pos or {obj.ox,obj.oy,obj.oz}
  local ox,oy,oz = unpack(pos)
  local xx,yy,zz = (ox+obj.x),(oy+obj.y),(oz+obj.z)
  local x,y,z = 0,0,0

  if string.find(tostring(L),"table:") then
    x,y,z = unpack(L)
  elseif L==nil then
    return ox,oy,oz,1
  else
    L="layer"..L
    x,y,z=gv(L..".x"),gv(L..".y"),gv(L..".z")
    r = radius or gv(L..".zoom")/2
    if gv(L..".zoom")==nil then
      return ox,oy,oz,1
      --x,y,z,r=0,0,0,0
    end
  end
  r = r or 200
  local xa,ya,za=xx-x,yy-y,zz-z
  local xb,yb,zb=x-xx,y-yy,z-zz

  if r<0 then
    xa,ya,za,xb,yb,zb = xb,yb,zb,xa,ya,za
  end
  r = math.abs(r)
  local l=(xa*xa + ya*ya + za*za)^.5
  local nx,ny,nz = xb/l,yb/l,zb/l
  local d = ( l<r and (1 - l/r) ) or 0
  d = d*f

if container==1 then
  if (l>r) then  --符号反転+feather=0だとContainer化
    ox = (xb-nx*r)*(f) + ox 
    oy = (yb-ny*r)*(f) + oy 
    oz = (zb-nz*r)*(f) + oz 
    if ret==nil or ret==0 then
      obj.ox,obj.oy,obj.oz = ox,oy,oz
    end
  end
  return ox,oy,oz,d
else
  if (l<r) then  --符号反転+feather=0だとContainer化
    ox = (xb-nx*r)*(f) + ox 
    oy = (yb-ny*r)*(f) + oy 
    oz = (zb-nz*r)*(f) + oz 
    if ret==nil or ret==0 then
      obj.ox,obj.oy,obj.oz = ox,oy,oz
    end
  end
  return ox,oy,oz,d
end

end

return Field

