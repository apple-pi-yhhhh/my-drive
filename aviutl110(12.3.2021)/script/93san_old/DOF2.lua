--[[
なんちゃってDepth of Field 偽被写界深度2 の関数版です。
exedit.aufと同じ場所にDOF2.luaを置くとスクリプト制御などで呼び出せます。
一度何らかの変数に関数を入れてから
local dof = require("DOF2")
実行してください。
dof(pos,blur,focusmode,focalpoint,aperture,alpha,fade,lens,mode)

直接使用することもできます。
require("DOF2")(nil,20)

すでにほかのオブジェクトで偽被写界深度2を使用してる場合、第二引数をnilにしておくと
偽被写界深度2のパラメーターと連動します。なので偽被写界深度2を実行ている場合はスクリプト制御に
require("DOF2")()
で動作します。


DOF2( [pos, blur, focusmode, focalpoint, aperture, fade , lens, mode])

	pos = {x,y,z} (エフェクトを掛ける対象の座標) 内部でオブジェクトの座標+相対座標が自動で加算されるので通常は nil か {0,0,0}でOKです
	blur = レンズブラーの強度、 偽被写界深度2を他のオブジェクトで使用している場合, nil にするとそちらと連動します。(0～100)程度
	focusmode = true だと焦点を目標に固定 (0はfalse扱いされます)
	focalpoint = 焦点距離     (-5000～5000)程度
	aperture = 焦点の合う範囲 (0～5000)程度
	fade = 第一戻り値の値(0～1の範囲) が変化する距離を％で指定します。100だとブラーと同様の変化です。 100より小さいと急激に減衰します。 (1～500)
	lens = 0だと通常のぼかしを使用、 デフォルトは1でレンズブラーです。
	mode = 1だと単純なカメラとの距離計算になります。 デフォルトは0

	第1戻り値に 透明度用の1~0の値, 第2戻り値にDepth値 0~ を返します。blurを0にすると戻り値だけ取り出せます。
	local alpha,depth = require("DOF2")(nil,0)


--使用例 (スクリーンにドットを描画するスクリプト)---------------------------


--縦横数
local nx,ny=16,9

--DOF2を呼び出し
local Dof = require("DOF2")

--ここでは元画像を仮想バッファに退避してループ毎にcopybufferで呼び出す仕込み
--画像を加工するエフェクトの場合重ね掛けにならない様にループ毎にloadかcopybufferで画像を初期化します。
if (obj.w+obj.h)/2 < 20 then
obj.copybuffer("tmp","obj")
else
obj.setfont("koruri",20)
obj.load("画像サイズが大きすぎます")
return
end

--カメラを向く
obj.setoption("billboard",3)

local sw,sh=obj.screen_w *1.0,obj.screen_h*1.0
for i=0,nx-1,1/2 do
for j=0,ny-1,1/2 do
--local frm=(i*ny*nz)*2+(j*nz)*2+k*2
local frm=(i*ny)*2+(j)*2
local R=rand(-1000,1000,-1,frm) *.001
local x=sw/(nx-1)*i - sw/2
local y=sh/(ny-1)*j - sh/2
local z=R*(obj.time)*50

-- 仮想バッファからオリジナル画像を呼ぶ
obj.copybuffer("obj","tmp")
-- エフェクトをかけつつ戻り値で透明度用の値を取得
local alp = Dof({x,y,z},10)
-- drawに透明度を適用
obj.draw(x,y,z, (R+1)/2 , alp)
end
end

------------------------------------

]]











local DOF2 = function(...)
 local T={...}
 local blur,focusmode,focalpoint,aperture,alpha,fade,lens,mode
 local pos = T[1] or {0,0,0} --{obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz}
 local objectpos = {obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz}

if DOF_TBL and T[2]==nil then
 	  blur,focusmode,focalpoint,aperture,alpha,fade,lens,mode = unpack(DOF_TBL)
		alpha = 0
else
	 blur       = T[2] or 0
	 focusmode  = T[3] or false
	 focalpoint = T[4] or 0
	 aperture   = T[5] or 0
	 alpha      = 0 --T[6] or 0
	 fade       = T[6] or 200
	 lens       = T[7] or 1
	 mode       = T[8] or 0
  if focusmode==0 then focusmode=false end
end


  local maxblur = blur * 2

	local Add = function(a,b)
		return {a[1]+b[1], a[2]+b[2], a[3]+b[3]}
	end

	local Sub = function(a,b)
		return {a[1]-b[1], a[2]-b[2], a[3]-b[3]}
	end

	local Scale = function(v,scale)
		return {v[1]*scale, v[2]*scale, v[3]*scale}
	end

	local Length = function(v)
		local x,y,z = v[1],v[2],v[3]
		return math.sqrt((x*x)+(y*y)+(z*z))
	end

	local Normalize = function(v)
		local x,y,z = v[1],v[2],v[3]
		local l = 1/math.sqrt(x*x + y*y + z*z)
		return {x*l, y*l, z*l}
	end

	local Dot = function(a,b)
		return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
	end

	local Cross = function(v0,v1)
		local x0,y0,z0 = v0[1],v0[2],v0[3]
		local x1,y1,z1 = v1[1],v1[2],v1[3]
		return {
		(y0 * z1) - (z0 * y1),
		(z0 * x1) - (x0 * z1),
		(x0 * y1) - (y0 * x1)
		}
	end

	local zoom = obj.getvalue("zoom")*.01
	local P0 = Add(pos,objectpos)
	local cam = obj.getoption("camera_param")
	local camT = {cam.tx,cam.ty,cam.tz}
	local camP = {cam.x,cam.y,cam.z}
	local camV = Sub(camT,camP)
	local screen = 1024 -- or cam.d

	--<< 中心点を使った回転を含みたくない場合はここを無効に//////////////////
	local ROT = function(p)
		local toRad=math.pi/180
		local x,y,z = p[1],p[2],p[3]
		local rx,ry,rz = obj.rx*toRad,obj.ry*toRad,obj.rz*toRad
		local x0=x*math.cos(rz)-y*math.sin(rz)
		local y0=x*math.sin(rz)+y*math.cos(rz)
		local z0=z*math.cos(ry)-x0*math.sin(ry)
		return {
			z*math.sin(ry)+x0*math.cos(ry),
			y0*math.cos(rx)-z0*math.sin(rx),
			y0*math.sin(rx)+z0*math.cos(rx)
		}
	end

	local C1 = ROT({obj.cx*zoom,obj.cy*zoom,obj.cz})
	      P0 = Sub(P0,C1)
	--///////////////////////////////////////////////////////ここまで >>

	local P1 = Sub(P0,camP)
	local d = Length(camV)
	if d<1 then d=1 end

	local L
	local N = Scale(camV,1/d) --Normalize(camV)
	local P = Scale(N,(focusmode and d) or screen)

	if mode==1 then
	 L = Length(Sub(camP,P0)) - (focusmode and d or screen )
	else
	 L = Dot(Sub(P1,P),N)
	end

  	local R = L - focalpoint
  	local D = math.abs(R)
	      D = (D-aperture)/screen
	if D<0 then D=0 end


	local fdp = fade*.01
	local alp = math.max(0,math.min(fdp-D,fdp)) / fdp
	if alp>1 then alp=1 end
	local Blur = math.abs(blur) * math.min(D,3)

	if (obj.alpha*obj.getvalue("alpha") > .0001 ) then

	if lens==0 then
	obj.effect("ぼかし","範囲",Blur)
	else
	obj.effect("レンズブラー","範囲",Blur,"サイズ固定",0)
	end

	end
	return alp,D,R
end

return DOF2
