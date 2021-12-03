
--[[

	スクリプト制御に一行書くだけでたくさん表示できるようにしたい用残念スクリプトいろいろ。
	CHAP_93.luaをexedit.aufと同じフォルダに配置して下さい。
	

	rikkyさんの拡張機能であるrikky_moduleが必要なのでお持ちでない方は
	http://hazumurhythm.com/wev/downloads/
	上記サイトでDLLパック(rikky_module)ver1.0 をダウンロードして
	zipファイルにあるrikky_module.dllをexedit.aufと同じフォルダに配置して下さい。
	
	また一部でeasing.luaを使用しているのでお持ちでない方はundufishさんの動画
	http://www.nicovideo.jp/watch/sm20813281
	上記動画urlからeasing_tra.zipをダウンロード、
	zipファイルにあるeasing.luaをexedit.aufと同じフォルダに配置して下さい。
	すでに持っている方もeasing.luaのみexedit.aufと同じフォルダに移動しておくと正しく動作します。


	最終更新2016/10/25

		-- 使い方(使わないと思いますが)

		-- スクリプト制御でluaを呼び出す
		require("CHAP_93")

		-- Rotcは回転行列を行います(回転がグループ制御のような挙動になる)
		-- 引数なしだと中心点を含んだオブジェクトの座標を返します。
		obj.ox, obj.oy, obj.oz = Rotc()
]]





require("rikky_module")
local type = rikky_module.type

--[[
local type = function(v)
  local s = tostring(v)
  if(s == v)then return "string" end
  if(s == "nil")then return "nil" end
  if(s == "true" or s == "false")then return "boolean" end
  if(string.find(s, "table:"))then return "table" end
  if(string.find(s, "function:"))then return "function" end
  if(string.find(s, "userdata:"))then return "userdata" end
  return "number"
end
]]



--単純なVector----------------------------------------------------------------------

local Sub	= function(a,b) return {a[1]-b[1],a[2]-b[2],a[3]-b[3]} end
local Add	= function(a,b) return {a[1]+b[1],a[2]+b[2],a[3]+b[3]} end
local Scale	= function(v,s) return {v[1]*s,v[2]*s,v[3]*s} end
local Mul	= function(a,b) return {a[1]*b[1],a[2]*b[2],a[3]*b[3]} end
local Div	= function(a,b) return {a[1]/b[1], a[2]/b[2], a[3]/b[3]} end
local Len	= function(v) return (v[1]*v[1]+v[2]*v[2]+v[3]*v[3]) end
local Lensq	= function(a,b) return Len( Sub(a,b) ) end
local Length	= function(v) return (v[1]*v[1]+v[2]*v[2]+v[3]*v[3])^.5 end
local Normalize	= function(v) local L=Length(v);local l=1/L; return {v[1]*l,v[2]*l,v[3]*l},L end
local Cross	= function(a,b) return  {a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1]} end
local Dot	= function(a,b) return (a[1]*b[1] + a[2]*b[2] + a[3]*b[3]) end
local eps	= 1e-08

-- 点Aと平面上の最近点(A=座標 P=平面上の点 N=平面の法線 )----------------------------
local NearPosOnPlane = function(A,P,N)
    local PA = {A[1]-P[1],A[2]-P[2], A[3]-P[3]}
    local d = (N[1]*PA[1] + N[2]*PA[2] + N[3]*PA[3])
    return {A[1]-(N[1]*d), A[2]-(N[2]*d), A[3]-(N[3]*d)}
end

-- 点Aと面abc上の最近点------------------------------------
local NearPosOnPoly=function(A,a,b,c)
  local AB = Sub(b,a)
  local BC = Sub(c,b)
  local N = Normalize(Cross(AB,BC))
  --local PA = Sub(A,a)
  --local d = Dot( N, PA )
  --return Sub(A,Scale(N,d)),N
  return NearPosOnPlane(A,a,N)
end

-- 点Pと線abの距離------------------------------------
local Dist_p_ab = function(P,a,b)
  local ab,aP = Sub(b,a), Sub(P,a)
  local cross   = Cross(ab,aP)
  local len   = Length(ab)
  local dst   = Length(cross) / len
  return dst
end

-- 点Pと線分abの距離------------------------------------
local Dist_p_isab = function(P,a,b)
	local AB = Sub(b,a)
  if	 ( Dot(AB,Sub(P,a) ) < eps ) then
    return Length(Sub(a,P))
  elseif ( Dot(AB,Sub(P,b) ) > eps ) then
    return Length(Sub(b,P))
  else
    return Dist_p_ab(P,a,b)
  end
end


-- 点Pと線ab上の最近点------------------------------------
local Pos_p_ab=function(P,a,b)
	local AB = Sub(b,a)
	local N = Normalize(AB) 	--線abの単位ベクトル
	local D = Dot(N, Sub(P,a))	--aPベクトルと内積
	return Add(a,Scale(N,D))
end

-- 点Pと線分ab上の最近点------------------------------------
local Pos_p_isab=function(P,a,b)
	local AB = Sub(b,a)
  	if	( Dot(AB,Sub(P,a) ) < eps ) then
   		return a
 	elseif	( Dot(AB,Sub(P,b) ) > eps ) then
   		return b
 	else
		return Pos_p_ab(P,a,b)
	end
end


-- 線分ABと平面の交点------------------------------------
local Intersect_plane_Line = function(
	A,-- 線分始点
	B,-- 線分終点
	n,-- 平面法線
	d,-- 平面法線の長さ
	PL-- = {n[1],n[2],n[3],d=d} -- ax+by+cz-d=0
  	)
	--平面上の点P
  	local P = Scale(n,d)
	local PA = Sub(A,P)
	local PB = Sub(B,P)
	--平面法線と内積
	local dot_PA = Dot(PA,n)
	local dot_PB = Dot(PB,n)
	--線端が平面上にあった時の誤差を0に
	if math.abs(dot_PA) < eps then  dot_PA = 0 end
	if math.abs(dot_PB) < eps then  dot_PB = 0 end
	--交差判定
	if (dot_PA == 0) and (dot_PB == 0) then
	-- 線端が平面上で計算不可
		return false
	elseif  ((dot_PA >= 0) and (dot_PB <= 0)) or ((dot_PA <= 0) and (dot_PB >= 0))  then
	-- 内積正負が異なれば交差
		local AB = Sub(B,A)
	-- 交点とAの距離 交点とBの距離 = dot_PA , dot_PB
		local ratio = math.abs(dot_PA) / ( math.abs(dot_PA) + math.abs(dot_PB) )
		return {
			A[1] + ( AB[1] * ratio ),
			A[2] + ( AB[2] * ratio ),
			A[3] + ( AB[3] * ratio )
		}
	else
	--交点なし
		return false
	end

end

--線AB,線CDで構成される２直線の交点(なければ最近点)------------------------------------
local PointOfIntersection2line = function(A,B,C,D)
  local AB = Sub(A,B)
  local CD = Sub(C,D)
  --計算不可
  if( Len(AB)==0) or (Len(CD)==0) then return 0,nil,nil end

  local n1 = Normalize(AB)
  local n2 = Normalize(CD)
  local w1 = Dot( n1, n2 )
  local w2 = 1 - w1*w1
  if( w2 == 0 ) then  return 0,false,false end
  local AC = Sub(A,C)
  local d1 = (Dot(AC,n1)-w1*Dot(AC,n2)) / w2
  local d2 = (w1*Dot(AC,n1)-Dot(AC,n2)) / w2
  local ret1,ret2
  --AB上の最近点
  ret1 = Add(A,Scale(n1,d1))
  --BC上の最近点
  ret2 = Add(C,Scale(n2,d2))

  if( Len(Sub(ret1,ret2)) < eps ) then
      return 1,ret1,ret2 --交点
  else
      return 2,ret1,ret2 --交点なし、最近点
  end
end


--法線との反射ベクトル s=法線(正規化),v=入射ベクトル------------------------------------
local Refrect_N=function(s,v)
  local t = -(s[1]*v[1] + s[2]*v[2] + s[3]*v[3])/(s[1]*s[1] + s[2]*s[2] + s[3]*s[3])
  return {v[1]+(t*s[1]*2), v[2]+(t*s[2]*2), v[3]+(t*s[3]*2)}
end



Vector = {
	Sub	 = Sub,
	Add	 = Add,
	Scale	 = Scale,
	Len	 = Len,
	Lensq	 = Lensq,
	Length	 = Length,
	Norm	 = Normalize,
	Normalize= Normalize,
	Cross	 = Cross,
	Dot	 = Dot,
	Mul	 = Mul,
	Div	 = Div,
	Pos_p_ab 	 = Pos_p_ab,
	Pos_p_isab 	 = Pos_p_isab,
	Refrect_N	 = Refrect_N,
	NearPosOnPlane 	 = NearPosOnPlane,
	NearPosOnPoly	 = NearPosOnPolye,
	Dist_p_ab	 = Dist_p_ab,
	Dist_p_isab 	 = Dist_p_isab,
	Intersect_plane_Line 	 = Intersect_plane_Line,
	PointOfIntersection2line = PointOfIntersection2line
	

}


-------------------------------------------------------------------
--[[
local Randmv = function(columns,tHold,tMin,tMax,w,h)

	--columns = 15;	 --number of columns in grid
	--tHold  = 1;	 --hold time (must be less than tmin)
	--tMin = 1;	 --minimum cycle time (can't be zero)
	--tMax  = 3;	 --maximum cycle time

columns,tHold,tMin,tMax = columns or 15,tHold or 1,tMin or 1,tMax or 3
w,h = w or obj.screen_w, h or obj.screen_h
local gap = w/columns;
local origin = {gap,gap,0};
local xGrid = columns - 1;
local yGrid = math.floor(h/gap) - 1;

local start = 0;
local ending = 0;
local j = 1;

while (time >= ending) do
  j = j+1;
  start = ending;
  ending = ending + rand(tMin,tMax,j,1);
end

local targetX = math.floor(rand(0,xGrid,j,1));
local targetY = math.floor(rand(0,yGrid,j,1));
local seed = (j-1);
local x = rand(0,100,j,1)*.01 ; --this is a throw-away value
local oldX = math.floor(rand(0,xGrid,seed,1));
local oldY = math.floor(rand(0,yGrid,seed,1));

	if(targetX == oldX and targetY == oldY) then
		return Add(origin, {oldX*gap,oldY*gap,0})
	elseif (time - start < tHold) then
		 return Add(origin, {oldX*gap,oldY*gap,0})
	else
		deltaX = math.abs(targetX - oldX);
		deltaY = math.abs(targetY - oldY);
 		xTime = (ending - start - tHold)*(deltaX/(deltaX + deltaY));
 		yTime = (ending - start - tHold)*(deltaY/(deltaX + deltaY));
 		 if (time < start + tHold + xTime) then
  	 		startPos = Add(origin, {oldX*gap,oldY*gap,0})
   	 		targetPos = Add(origin, {targetX*gap,oldY*gap,0})
   	 		return {
				require("easing").inQuad((time - start - tHold)/xTime, startPos[1], targetPos[1], targetPos[1]-startPos[1]),
				require("easing").inQuad((time - start - tHold)/xTime, startPos[2], targetPos[2], targetPos[2]-startPos[2])
			}
 		else
  	 		startPos  = Add(origin , {targetX*gap, oldY*gap,0}})
   			targetPos = Add(origin , {targetX*gap, targetY*gap,0})
   			 return {
				require("easing").outQuad((time - start - tHold - xTime)/yTime, startPos[1], targetPos[1],targetPos[1]-startPos[1]),
				require("easing").outQuad((time - start - tHold - xTime)/yTime, startPos[2], targetPos[2],targetPos[2]-startPos[2])
			}
  		end
	end
end

]]

--obj.loadで初期化されるデータのコピー、ロード-----------------------------------------------------
	--	local each_info = Param()	-- 引数なしでパラメータを返す
	--	obj.load(...)			-- obj.loadを使うとパラメータは初期化されてしまうので
	--	Param(each_info) 		-- 元のパラメータをロード
Param = function(param)
	if tostring(param):find("table:") then
		obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.cx,obj.cy,obj.cz,obj.aspect=unpack(param)
	elseif param==0 then
		obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.aspect,obj.cx,obj.cy,obj.cz = 0,0,0,1,1,0,0,0,0
	else
		return {obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.cx,obj.cy,obj.cz,obj.aspect}
	end
end


--他のレイヤー座標をまとめる---------------------------------------------------

GL = function(...)
	local tx = {[0]=".x",".y",".z"}
	local V,A,n = {},{},3
	for k=1,select("#",...) do
		A[k]={}
		for i=0,n-1 do
			local val = obj.getvalue("layer"..select(k,...)..tx[i])
			V[k*n+i-n+1] = val
			A[k][i+1] = val
		end
	end
	return V,A
	
end

GL2 = function(...)
	local tx = {[0]=".x",".y",".z"}
	local V,A,n = {},{},3
	for k=1,select("#",...) do
		A[k]={}
		for i=0,n-1 do
			local val = obj.getvalue("layer"..select(k,...)..tx[i])
			--V[k*n+i-n+1] = val
			A[k][i+1] = val
		end
	end
	return A
end

--tが1増える毎にテーブルvの値を線形補完(加減速)で巡回-----------------------------------------------
Swap_value = function(t,v)
	local v1=v
	local l=#v
	local i,d = math.modf(t)
	d = (math.cos(math.pi*d)-1)*-.5
	local x,y,z
	x = v[i%l+1]*(1-d) + v1[(i+1)%l+1]*d
	return x
 end

--tが0から1に増えるとカラーコードcola,colbを交換-----------------------------------------------
Swap_col = function(t,cola,colb)
	local r,g,b=RGB(cola)
	local R,G,B=RGB(colb)
	return RGB(r*(1-t)+R*t, g*(1-t)+G*t, b*(1-t)+B*t)
end

--テーブルaとbの値を入れ替え----------------------------------------------------------------------
--tが0から1に増えると、 aとbを線形補完で交換
--a,b はどちらがか数値、テーブルでも可。
--戻り値のタイプ、テーブルの長さはaに従う。(ただしt==0ならaを、t==1ならbをそのまま返す)

Swap_table = function(t,a,b)
	if (t==0) then
		return a
	elseif (t==1) then
		return b
	else
		t = (math.cos(math.pi*t)-1)*-.5
		local flag=0
 
      		if not tostring(a):find("table:") then
      			a  ={a}
      			flag = 1
      		end

      		if not tostring(b):find("table:") then
      			b = {b}
      		end

        	if (flag==1) then
        		return (a[1]*(1-t) + b[1]*t)
       		else
          		local c = a
           		for kb,vb in pairs(b) do
            		for ka,va in pairs(a) do
              			if (kb==ka) then
                			c[kb] = va*(1-t) + vb*t
              			end
           		end
           		end
          		return c
		end
        end

end

--テーブルの結合----------------------------------------------------------------------

Merge_table = function(...)
	local t = {...}
	local T = t[1]
	for i=2,#t do
	for j=1,#t[i] do
		table.insert(T,t[i][j])
	end
	end
	return T
end


--linear補完----------------------------------------------------------------------
Linear = function(t, t1, t2, v1, v2, nolimit )
     v1 = v1 or 0
     v2 = v2 or 1
     local c = (t2 - t1)
     local n = t/c - t1/c
     local V = v2-v1
           V = V * n + v1
      if (nolimit==1) then
         return V
      else
       if (v1>v2) then v1,v2=v2,v1 end
         V = math.max(v1,math.min(v2,V))
         return V
      end
  end


--progressが0～1の動作で0～1になるずれた値が入ったテーブルを返します。----------------------------------------------------------------
 
--順番の入れ替えは戻り値をShuffle関数へ。

--[[
Progress = function(progress, overlap, num, ease)
	overlap = overlap or 1
	num = num or 10
	local total = 1 + overlap*(num-1)
	local t = progress * total
	local T = {}
	for i=0,num-1 do
		local v = t - overlap*i 
		v = (v<0 and 0) or (v>1 and 1) or v
			if ease and #ease>0 then
				v = require("easing")[ease](v,0,1,1)
			end
		T[i+1] = v
	end
	return T
end
]]

Progress = function(progress, wait, num, ease)
	wait = wait or 1
	num = num or 10
	local total_dulation = 1 + wait*(num-1)
	local t = progress * total_dulation
	local T = {}
	-- easeing --------------------------------------
	if ease and #tostring(ease)>0 then
		local ease_s = tostring(ease)
		local E = require("easing")
		if (ease_s):find("%d") then
			if tonumber(ease_s)<=41 then 
				local ez = {
				"linear",
				"inSine","outSine","inOutSine","outInSine",
				"inQuad","outQuad","inOutQuad","outInQuad",
				"inCubic","outCubic","inOutCubic","outInCubic",
				"inQuart","outQuart","inOutQuart","outInQuart",
				"inQuint","outQuint","inOutQuint","outInQuint",
				"inExpo","outExpo","inOutExpo","outInExpo",
				"inCirc","outCirc","inOutCirc","outInCirc",
				"inElastic","outElastic","inOutElastic","outInElastic",
				"inBack","outBack","inOutBack","outInBack",
				"inBounce","outBounce","inOutBounce","outInBounce"
				}
				ease = ez[tonumber(ease_s)] 
			end
		end

		for i=0,num-1 do
			local v = t - wait*i 
			v = (v<0 and 0) or (v>1 and 1) or v
			v = E[ease](v,0,1,1)
			T[i+1] = v
		end
	else
	--------------------------------------------------
		for i=0,num-1 do
			local v = t - wait*i 
			v = (v<0 and 0) or (v>1 and 1) or v
			T[i+1] = v
		end
	end

	return T
end

--<< Shuffle (num,order,seed) >>---------------------------------------------------------

Shuffle = function(num,order,seed,scl)

--[[
orderはさつきさんのTAにある順と同じです。

numが数値の場合はシャッフルされた値が入ったテーブルをnum要素分作って返します。
戻り値のテーブルをfor分などで使用
	 num = 12
	 order = 0
	 t = Shuffle(num,order,0)
	 {1,2,3,4,5,6,7,8,9,10,11,12} --長さがnumのテーブル

	 num = 5
	 order = 2 -- order2はランダム
	 t = Shuffle(num,order,0)
	 {5,1,3,4,2}

numがテーブルだった場合はテーブル自体をシャッフルして返します。
第4引数は値をスケールします。(テーブルの中身が数値の場合)
	num = {1,2,3,4,"A"}
	order = 2
	t = Shuffle(num,order,0)
	 {"A",1,3,4,2}

orderが3，4，5の場合は要素が抜けます(半分で折り返すため)
	num = {1,2,3,"A",5,6,7}
	order = 5
	t = Shuffle(num,order,0)
	{1,2,3,"A",3,2,1}

]]
local tbl
if tostring(num):find("table:") then tbl=num; num=#tbl  end
order=order or 0
seed=seed or -1
local index={}
  for i=0,num-1 do
    local k=i+1
    if(order<1) then
      index[k] =i 
    elseif(order<2) then
      index[k]=num-1-i 
    elseif(order<3) then
      local es={}
      for j=0,num-1 do
        es[j+1]=j
      end
      for j=0,num-1 do
        local dest = 0
        dest=rand(0,num-1, -num - math.abs(seed),j+1)
        local swap=es[j+1]
        es[j+1]=es[dest+1]
        es[dest+1]=swap
      end
      index[k]=es[k]
    elseif(order<4) then
	index[k]=math.floor(rand(0,100*(num-1),seed,i)*.01 +.5) 
    elseif(order<5) then
	index[k]=math.floor(math.abs((num-1)/2-i) )*2
    else
        index[k]=( (num-1)/2-math.abs((num-1)/2-i) )*2
    end
  end

if tbl then
  local t=tbl
  tbl={}
  if scl then
 	 for i=1,num do
   	 	local j=index[i]+1
   	 	tbl[i] = t[j]*scl
 	 end
  else
  	for i=1,num do
    		local j=index[i]+1
   		tbl[i] = t[j]
 	end
  end
  return tbl
else
return index
end
end

--回転行列----------------------------------------------------------------------

Rot = function(v,r)
	v = v or {0,0,-1}
	r = r or {obj.rx,obj.ry,obj.rz}
	local tR = math.pi/180
	local cos,sin=math.cos,math.sin
	local rx,ry,rz = r[1]*tR, r[2]*tR, r[3]*tR
	local x,y,z = v[1],v[2],v[3]
	local x0=x*cos(rz)-y*sin(rz)
	local y0=x*sin(rz)+y*cos(rz)
	local z0=z*cos(ry)-x0*sin(ry)
	return z*sin(ry)+x0*cos(ry), y0*cos(rx)-z0*sin(rx), y0*sin(rx)+z0*cos(rx)
end

--相対座標、中心点を移動後に回転行列----------------------------------------------------------------------
Rotc = function(pos,anc,rot)
	pos = pos or {obj.ox, obj.oy, obj.oz}
	anc = anc or {obj.cx, obj.cy, obj.cz}
	rot = rot or {obj.rx, obj.ry, obj.rz}
	local zoom = obj.getvalue("zoom")*.01
	local ox,oy,oz = pos[1],pos[2],pos[3]
	local cx,cy,cz = anc[1]*zoom,anc[2]*zoom,anc[3]*zoom
	local x,y,z = Rot({ox-cx, oy-cy, oz-cz},rot)
	return  x+cx, y+cy, z+cz
end

--補完付きランダム----------------------------------------------------------------------
Shake = function(interval,min,max,seed,frm)
      if not min then min,max,seed=-1,1,0 end
        frm = frm or obj.frame
        local min,max=min*1000,max*1000
        local t = obj.time*100/interval
        local p = math.floor(t)
        t = t-p
        local pl={}
        pl[1]=t
        for i=0,3*4-1 do
          local s = seed
          pl[i+2]=rand(min,max,s,p+math.floor(i/4)+frm)*.001
        end
        return obj.interpolation(unpack(pl))
end



--カメラ深度エフェクト----------------------------------------------------------------------
local Depthfx = function(
	pos,		-- オブジェクト座標{x,y,z}
	focalpoint,	-- 焦点の前後
	startfade,	-- フェード開始位置
	vanish,		-- フェード開始からの消滅距離
	near_startfade, -- フェード開始位置(焦点より手前)
	near_vanish,	-- フェード開始からの消滅距離 (焦点より手前)
	focusmode	-- 1で焦点を目標点に固定
	)
	pos = pos or {obj.x+obj.ox,obj.y+obj.oy,obj.z+obj.oz}
	focalpoint,startfade,vanish = focalpoint or 0,startfade or 100,vanish or 2500
	near_startfade,near_vanish = near_startfade or startfade,near_vanish or vanish
	local mul=function(v,s) return {v[1]*s,v[2]*s,v[3]*s} end
	local dot=function(a,b) return (a[1]*b[1]+a[2]*b[2]+a[3]*b[3]) end
	local c = obj.getoption("camera_param")
	local e = {c.tx-c.x, c.ty-c.y, c.tz-c.z}
	local d = (e[1]*e[1]+e[2]*e[2]+e[3]*e[3])^.5
	local n = mul(e,1/d)
	local pd = (focusmode) and (d+focalpoint) or (1024+focalpoint)
	local pl = mul(n,pd)
	local pv = {pos[1]-c.x, pos[2]-c.y, pos[3]-c.z}
	local depth = dot({pv[1]-pl[1],pv[2]-pl[2],pv[3]-pl[3]},n)
	if depth<0 then startfade,vanish = near_startfade,near_vanish end
	local D = math.abs(depth)-startfade
	D = D<0 and 0 or D
	D = D>vanish and vanish or D
	D = (1-D/vanish)
	return D,depth,n,pd
end

--長方形作成----------------------------------------------------------------------
Rect = function(w,h,col,alp,x,y,rot,mode)
	local Rz = function(x,y,r)
		r = r*math.pi/180
		return x*math.cos(r)-y*math.sin(r), x*math.sin(r)+y*math.cos(r)
	end
	x,y,rot=x or 0,y or 0,rot or 0
	rot=rot or 0
	x=x+.5
	y=y+.5
	w,h=w/2,h/2
	local _w,_h=-w,-h
	if mode==1 then w,h,_w,_h=w*2,h*2,0,0 end
	if col then obj.putpixel(0,0,col,1) end
	local x0,y0=Rz(_w+x,_h+y,rot)
	local x1,y1=Rz( w+x,_h+y,rot)
	local x2,y2=Rz( w+x, h+y,rot)
	local x3,y3=Rz(_w+x, h+y,rot)
	obj.drawpoly(
		x0,y0,0,
		x1,y1,0,
		x2,y2,0,
		x3,y3,0,
		0,0, 0,0, 0,0, 0,0, alp or 1
	)
	
	return {{x0,y0,0},
		{x1,y1,0},
		{x2,y2,0},
		{x3,y3,0}
	}
end

Rect2 = function(w,h,col,alp,x,y,rot)
	Rect(w,h,col,alp,x,y,rot,1)
end

--台形作成----------------------------------------------------------------------
--上下の辺,左右の辺の長さを入力して台形を作成します。
-- w,hはテーブルで w={上辺,底辺} h={左辺,右辺}と別に指定できます。
-- 長方形や三角形もこれで代用可能
Trapez = function(w,h,col,alp,x,y,rot,shear,scl)
	x,y,rot,shear,scl= x or 0, y or 0,rot or 0, shear or 0,scl or 1
	if not tostring(w):find("table:") then
		w = {w,w}
	end
	if not tostring(h):find("table:") then
		h = {h,h}
	end
	tw  = w[1]/2*scl
	bw  = w[2]/2*scl
	lh  = h[1]/2*scl
	rh  = h[2]/2*scl
	local Rz = function(x,y,r)
		r = r*math.pi/180
		return x*math.cos(r)-y*math.sin(r), x*math.sin(r)+y*math.cos(r)
	end
	if col then obj.putpixel(0,0,col,1) end
	local x0,y0 = Rz(-tw + x+shear,-lh+y,rot)
	local x1,y1 = Rz( tw + x+shear,-rh+y,rot)
	local x2,y2 = Rz( bw + x-shear, rh+y,rot)
	local x3,y3 = Rz(-bw + x-shear, lh+y,rot)
	obj.drawpoly(x0,y0,0, x1,y1,0, x2,y2,0, x3,y3,0, 0,0, 0,0, 0,0, 0,0, alp or 1)
	return { {x0,y0,0}, {x1,y1,0}, {x2,y2,0}, {x3,y3,0} }
end

Drawline = {}
--点p0と点p1を結ぶ直線(2D)----------------------------------------------------------------------

P_line = function(
	p0,	--座標 {x,y}
	p1,	--座標 {x,y}
	width,	--[線幅]
	col,	--[色]
	alp,	--[透明度]
	st,	--[消滅開始距離]
	va,	--[消滅までのフェード範囲]
	t	--[0~1で線を伸ばす]
	)
	width = width or 1
	width = width*.5
	st = st or 1500
	va = va or 2500
	local x0,y0=p0[1],p0[2]
	local x1,y1=p1[1],p1[2]
	local x,y=(x1-x0),(y1-y0)
	local L=(x*x+y*y)^.5

	if L>(st+va) then
		return 0
	else
		t = t or 1
		t = math.max(0,math.min(1,t))
		local mul=function(v,s) return {v[1]*s,v[2]*s} end
		local add=function(a,b) return {a[1]+b[1],a[2]+b[2]} end
		p1 = add( mul(p1,t), mul(p0,1-t))
		x1,y1 = p1[1],p1[2]
  		local l = L-st
 		l = l<0 and 0 or l
 		l = l>va and va or l
		l = (1-l/va)^2
 		width = width*l
 		local xc,yc= -(y/L)*width, (x/L)*width
  		if col then obj.putpixel(0,0,col,1) end
  		obj.drawpoly(
  			x0+xc,y0+yc,0,
   			x1+xc,y1+yc,0,
   			x1-xc,y1-yc,0,
   			x0-xc,y0-yc,0,
   			0,0, 0,0, 0,0, 0,0,alp or 1
 	 	)
		return 1
 	end
end

--p0,p1を結ぶ直線(3D,カメラ用)----------------------------------------------------------------------

P_line3D = function(
	p0,	--座標 {x,y,z}
	p1,	--座標 {x,y,z}
	width,	--[線幅]
	col,	--[色]
	alp,	--[透明度]
	st,	--[消滅開始距離]
	va,	--[消滅までのフェード範囲]
	nst,	--[消滅開始距離(st以下)]
	nva,	--[消滅までの範囲(nst以下～0まで)]
	t	--[0~1で線を伸ばす]
	)
	local V = Vector
	width = width or 1
	alp = alp or 1
	st,va = st or 500,va or 1000
	nst,nva = nst or 1 ,nva or 0
	t = t or 1
	t = math.max(0,math.min(1,t))
	p0[3] = p0[3] or 0
	p1[3] = p1[3] or 0
	local a=V.Sub(p1,p0)
	local len=V.Length(a)
	if len>(st+va) then return 0,p1 end
	if len<(nst-nva) then return 0,p1 end

	local lin = Linear(len,st,va,1,0)
	lin = lin*Linear(len,nva,nst ,0,1)
	width = width*lin
	alp = alp*lin
	local g=obj.getvalue

	if col then obj.putpixel(0,0,col,1) end

	local c=obj.getoption("camera_param")
	local b={c.x-p0[1], c.y-p0[2], c.z-p0[3]}
	local n = V.Cross(a,b)
	local l = V.Length(n)
	local nx,ny,nz = (n[1]/l)*width*.5, (n[2]/l)*width*.5 ,(n[3]/l)*width*.5
	p1 = V.Add( V.Scale(p1,t), V.Scale(p0,1-t))
	obj.drawpoly(
		p0[1]-nx,p0[2]-ny,p0[3]-nz,
		p1[1]-nx,p1[2]-ny,p1[3]-nz,
		p1[1]+nx,p1[2]+ny,p1[3]+nz,
		p0[1]+nx,p0[2]+ny,p0[3]+nz,
		0,0,0,0,0,0,0,0,alp
		)
	return 1
end

--座標の入ったテーブルtを受け取って線で結ぶ----------------------------------------------------------------------
-- mode = 0 :順 	maxnum :増えると順に引く、nilで一周。
-- mode = 1 :すべて巡回 maxnum :増えると順に引く、ネズミ算式に増える。テーブルで第二引数で点あたりの最大数を制限 ・ 省略時は1(mode 0と同じ)
-- mode = 2 :すべて巡回 maxnum :点一つ当たり何本引くか指定のみ。
-- tは {{x,y,z},{x,y,z},{x,y,z}...}の形式。 一応 {x,y,z,x,y,z,x,y,z...}でも可。

C_line = function(mode,t,max_num,w,col,alp,st,va)
	w,col,alp,st,va = w or 1, col or 0xffffff, alp or 1, st or 2000, va or 2500
	if (max_num==nil) then
		max_num = {#t,1}
	elseif not tostring(max_num):find("table:") then
		max_num = {max_num, 1}
	end
	local maxnum = max_num[1]
	local maxcount = max_num[2] or 1
	if maxcount<0 then maxcount=#t end
	maxnum,af = math.modf(maxnum+1)
	--maxnum = math.min(maxnum,#t+0)
	af = (math.cos(math.pi*af^2)-1)*-.5
	--require("rikky_module")
	--local r_c = rikky_module.camerainfo

	if not tostring(t[1]):find("table:") then
		local t0={}
		for i=1,#t/3 do
			t0[i]={t[i*3-2],t[i*3-1],t[i*3]}
		end
		t = t0
	end
	--local p1 = t[1] --第二戻り値(線の先端座標)
	if mode==0 then
		for i=1,#t do
			local progress = 1
			if maxnum<i then 
				progress = 0
			elseif maxnum==i then 
				progress = af
			end
			n = P_line3D(t[i],t[i%#t+1],w,col,alp,st,va,nil,nil,progress)
		end
	elseif (mode==1) then

		for i=1,#t-1 do
			local n = 0
			for j=i,#t do
				local progress = 1
				if maxnum<j then 
					progress = 0
				elseif maxnum==j then 
					progress = af
				end
				--local nc,p1 = P_line3D(t[i],t[j%#t+1],w,col,alp,st,va,nil,nil,progress) 
				n = n + P_line3D(t[i],t[j],w,col,alp,st,va,nil,nil,progress) 
				if (n>=maxcount) then break end
			end
		end
	elseif (mode==2) then
		-- local maxcount,af = math.modf(maxcount)
 		for i=1,#t-1 do
			local n = 0
			for j=i,#t do
				local progress = 1				
				--if maxcount<n then 
				--	progress = 0
				--elseif j==maxcount then 
				--	progress = af
				--end
				n = n + P_line3D(t[i],t[j%#t+1],w,col,alp,st,va,nil,nil,progress)
				if (n>=maxnum) then break end
			end
		end
	end
	return p1
end

Drawline = {
		P2 = P_line,
		P3 = P_line3D,
		C  = C_line 
	}

--円形配置----------------------------------------------------------------------
--戻り値に配置の座標テーブル。 テーブルは{{x,y,z},{x,y,z}...} 形式
Cir = function(
		rad, 	--半径
		num, 	--描画数
		rot, 	--[全体のz軸回転 (nilだと中心を向かない) ]
		comp, 	--[0～100 で描画数を変更]
		zoom, 	--[拡大率 ]
		alp, 	--[透明度]
		rnd, 	--[半径のランダム]
		itv 	--[ﾗﾝﾀﾞﾑの移動間隔]
		)
	--obj.setoption("focus_mode","fixed_size")
	local is_tbl_r = tostring(rad):find("table:")
	local is_tbl_a = tostring(alp):find("table:")
	if is_tbl_r then num = #rad end
	local t = {}
	zoom = zoom or 1
	comp = comp or 100
	comp = math.min(1,comp*.01)
	local r = math.pi*2/num
	local r_flg = (rot==nil)
	if r_flg  then rot_=0 else rot_ = rot -90 end
	rot = rot or 0
	rot = rot*math.pi/180 - math.pi/2
	local rz = 360/num
	num = math.floor(num*comp)

	for i = 0,num-1 do
		local alpha=1
		if is_tbl_a then 
			alpha = alp[i%(#alp)+1]
		else
			alpha = alp or 1
		end
		local rz_ = rz*i
		local rd = 0
		if not (r_flg) then
			rz_=rz*i + rot_
		end
		if (rnd) then
			itv = itv or 200
			local iv=itv*rand(10,100,1,i+1)
			rd = rnd*Shake(itv+iv,-100,100,itv,i)*.01
		end
		local radius
		if (is_tbl_r) then
			radius = rad[i%#rad + 1 ] + rd
		else
 			radius = rad + rd
		end
		local x,y = math.cos(r*i+rot)*radius, math.sin(r*i+rot)*radius
		t[i+1]={x,y,0}
		if alpha>0 and zoom>0 then
			obj.draw(x,y,0,zoom,alpha,0,0,rz_)
		end
	end
	table.insert(t,t[1])
	return t
end

--[[

--円形配置----------------------------------------------------------------------
--戻り値に配置の座標テーブル。 テーブルは{{x,y,z},{x,y,z}...} 形式
Cir2 = function(
		rad, 	--半径
		num, 	--描画数
		rot, 	--[全体のz軸回転 (nil以外だと個別に中心を向く) ]
		comp, 	--[0～100 で描画数を変更]
		zoom, 	--[拡大率]
		alp, 	--[透明度]
		rnd, 	--[半径のランダム]
		itv 	--[ﾗﾝﾀﾞﾑの移動間隔]
		)
	obj.setoption("focus_mode","fixed_size")
	local scl = 100
	local is_tbl  = tostring(rad):find("table:")
	local is_tbla = tostring(alp):find("table:")
	if is_tbl then 
		n=#rad
		scl = num
	end
	scl = scl*.01
	local tb = {}
	zoom = zoom or 1
	comp = comp or 100
	comp = math.min(1,comp*.01)
	local r = math.pi*2/n
	local r_flg = (rot==nil)
	if r_flg  then rot_=0 else rot_ = rot -90 end
	rot = rot or 0
	rot = rot*math.pi/180 - math.pi/2
	local rz = 360/n
	local num = math.floor(n*comp)

	for i = 0,num-1 do
		local alpha = 1
		if is_tbla then 
			alpha = alp[(i+1)%(#alp)+1]
		else
			alpha = alp or 1
		end

		local rd=0
		if not (r_flg) then
			obj.rz=rz*i + rot_
		end
		if (rnd) then
			itv = itv or 200
			local iv=itv*rand(10,100,1,i+1)
			rd = rnd*Shake(itv+iv,-100,100,itv,i)*.01
		end
		local radius
		if (is_tbl) then
			radius = rad[(i+1)%#rad + 1 ]*scl + rd
		else
 			radius = rad + rd
		end
		local x,y = math.cos(r*i+rot)*radius, math.sin(r*i+rot)*radius
		tb[i+1]={x,y,0}
		obj.draw(x,y,0,zoom,alpha)
	end
	obj.rz=orz
	table.insert(tb,tb[1])
	return tb
end
]]


--グリッド配置----------------------------------------------------------------------

--[[
	local num_xyz = {10,10,3}
	lcoal div_xyz = {100,100,400}
	Dots(num_xyz, div_xyz)
]]

Dots = function(
		n, 	-- 描画数	( n or {xn,yn} or {xn,yn,zn} )
		d, 	-- 描画間隔	( d or {xd,yd} or {xd,yd,zd} )
		x, 	-- 描画位置ズラシx
		y, 	-- 描画位置ズラシy
		zoom, 	-- 拡大率
		alpha, 	-- 透明度
		rot, 	-- z軸の回転(個別)
		spfield	-- Spherical Field Option = {layer, radius , strength, container}
		)
	alpha,zoom,rot=alpha or 1,zoom or 1,rot or 0
	if not tostring(n):find("table:") then n={n,1,1} end
	if not tostring(d):find("table:") then d={d,d,1} end
	local nx,ny,nz=n[1],n[2],n[3] or 1
	local dw,dh,dz=d[1],d[2],d[3] or 0
	x,y=x or 0,y or 0
	local Field = function(v) return v[1],v[2],v[3] end
	local pos, radius , strength, container = 0,0,0,0
	if (spfield) then
		Field = SPField
		pos, radius , strength, container = unpack(spfield)
	end
	local fx,fy,fz=(nx-1)*dw/2,(ny-1)*dh/2,(nz-1)*dz/2
	local t = {}
	local t2 = {}
	--local t3 = {}
	local count=0
	for i=0,nx-1 do
		t2[i+1]={}
		--t3[i+1]={}
		for j=0,ny-1 do
			t2[i+1][j+1]={}
			--t3[i+1][j+1]={}
			for k=0,nz-1 do
				count = count+1
				local X=i*dw-fx + x
				local Y=j*dh-fy + y
				local Z=k*dz-fz
				X,Y,Z = Field({X,Y,Z},pos, radius , strength, container)
				t[count] = {X,Y,Z}
				t2[i+1][j+1][k+1] = {X,Y,Z}
				local Xr,Yr,Zr = Rot({X,Y,Z})
				--t3[i+1][j+1][k+1] = {Xr,Yr,Zr}
				obj.draw(Xr,Yr,Zr,zoom,alpha,-obj.getvalue("rx"),-obj.getvalue("ry"),-obj.getvalue("ry")+rot)
			end
		end
	end
	t.v  = t2
	--t.vr = t3
	return t
end

--球状配置----------------------------------------------------------------------

Sphere = function(rad,num,zoom,alpha,scale,spfield)
	zoom,alpha = zoom or 1,alpha or 1
	scale = scale or {1,1,1}
  	if not tostring(num):find("table:") then
		num = {num,math.ceil(num/2)}
  	end
	if not tostring(scale):find("table:") then
		scale = {scale,scale,scale}
	end
	local Draw = function(...) return 0 end
	if alpha>0 and zoom>0 then Draw=obj.draw end
	local Field = function(v) return v[1],v[2],v[3] end
	local pos, radius , strength, container = 0,0,0,0
	if (spfield) then
		Field = SPField
		pos, radius , strength, container = unpack(spfield)
	end
	
	local wn ,hn = num[1],num[2]
	local sx ,sy, sz = scale[1] or 1,scale[2] or 1,scale[3] or 1
	local p = math.pi*2/wn
	local t = math.pi/hn
	local sin,cos,toRad = math.sin,math.cos,math.pi/180
	local T  = {}
	local T2 = {}
	local T3 = {}
	local cn = 0
	for i=0,wn-1 do
		T2[i+1]={}
		T3[i+1]={}
		for j=0,hn do
			cn = cn+1
			local x,y,z = sin(t*j)*cos(p*i)*rad *sx, -(cos(t*j)*rad + math.pi/2)*sy , sin(t*j)*sin(p*i)*rad * sz
			x,y,z = SPField({x,y,z},pos, radius , strength, container)
			T2[i+1][j+1]={x,y,z}
			T[cn] = {x,y,z}
			x,y,z = Rot({x,y,z})
			T3[i+1][j+1]={x,y,z}
			Draw(x,y,z,zoom,alpha,-obj.getvalue("rx"),-obj.getvalue("ry"),-obj.getvalue("rz"))
		end
	end
	T.v  = T2
	T.vr = T3
	return T
end

--球状配置2----------------------------------------------------------------------

Sphere2 = function(rad,num,zoom,alpha,scale,spfield,figure,figsize,deform)
	zoom,alpha = zoom or 1,alpha or 1
	scale = scale or {1,1,1}
  	if not tostring(num):find("table:") then
		num = {num,math.ceil(num/2)}
  	end
	if not tostring(scale):find("table:") then
		scale = {scale,scale,scale}
	end
	local Draw = function(...) return 0 end
	local Load = function() return end
	if alpha>0 and zoom>0 then Draw=obj.draw Load=obj.load end
	local Field = function(v) return v[1],v[2],v[3] end
	local pos, radius , strength, container = 0,0,0,0
	if (spfield) then
		Field = SPField
		pos, radius , strength, container = unpack(spfield)
	end
	
	local wn ,hn = num[1],num[2]
	local w,h = obj.getpixel()
	local dw,dh = (w-1)/wn,(h-1)/(hn+1)
	local C = {}
	deform = deform or 0
	obj.pixeloption("type","rgb")
	for i=0,wn-1 do
		C[i+1]={}
		for j=0,hn do
			local r,g,b,a = obj.getpixel(dw*i,dh*j)
			C[i+1][j+1]={r=r,g=g,b=b,alpha=a, Y=(0.29891 * r + 0.58661 * g + 0.11448 * b) /255, col=RGB(r,g,b)}
		end
	end
	local sx ,sy, sz = scale[1] or 1,scale[2] or 1,scale[3] or 1
	local p = math.pi*2/wn
	local t = math.pi/(hn)
	local sin,cos,toRad = math.sin,math.cos,math.pi/180
	local T  = {}
	local T2 = {}
	local T3 = {}
	local cn = 0
	obj.setoption("billboard",3)
	for i=0,wn-1 do
		T2[i+1]={}
		T3[i+1]={}
		for j=0,hn do
			local m = C[i+1][j+1]
			cn = cn+1
			local rad = rad + rad * m.Y * deform
			local x,y,z = sin(t*j)*cos(p*i)*rad *sx, -(cos(t*j)*rad + math.pi/2)*sy , sin(t*j)*sin(p*i)*rad * sz
			x,y,z = SPField({x,y,z},pos, radius , strength, container)
			T2[i+1][j+1]={x,y,z,pix=m}
			T[cn] = {x,y,z}
			x,y,z = Rot({x,y,z})
			T3[i+1][j+1]={x,y,z}
			
			local col = m.col
			obj.load("figure",figure or "円",col,figsize or 10)
			Draw(x,y,z,zoom ,alpha*m.alpha/255 ,-obj.getvalue("rx"),-obj.getvalue("ry"),-obj.getvalue("rz"))
		end
	end
	obj.pixeloption("type","col")
	T.v  = T2
	T.vr = T3
	return T
end

--トーラス配置----------------------------------------------------------------------
--戻り値に配置の座標テーブル。 戻り値のテーブルは  t[i] = {x,y,z} 形式 と t[i][j] = {x,y,z} (i=大半径,j=小半径)
Torus = function(NUM,num,R,r,P,T,ROT,rot,alpha,zoom,spfield,shakeparam)
	local Tor = function(R,r,p,t)
		local x = R*math.cos(t) + r*math.cos(p)*math.cos(t)
		local y = R*math.sin(t) + r*math.cos(p)*math.sin(t)
		local z = r*math.sin(p)
		return x,y,z
	end
	local Shakev= function(v) return v[1],v[2],v[3] end
	local Field = function(v) return v[1],v[2],v[3] end
	local pos, radius , strength, container = 0,0,0,0
	if (spfield) then
		Field = SPField
		pos, radius , strength, container = unpack(spfield)
	end
	if (shakeparam) and (#shakeparam>4) then
		Shakev = function(v,cn,shakeparam)
			local itv ,min,max,seed,frm,am = unpack(shakeparam)
			am=am or 100
			am = am*.01
			local x = v[1] + Shake(itv+1 ,min,max, seed+100,cn+frm+1000)*am
			local y = v[2] + Shake(itv+2 ,min,max, seed+200,cn+frm+2000)*am
			local z = v[3] + Shake(itv+3 ,min,max, seed+300,cn+frm+3000)*am
			return x,y,z
		
		end
	end
	
	local pd = math.pi/180
	local rx,ry,rz=obj.getvalue("rx"),obj.getvalue("ry"),obj.getvalue("rz")
	alpha,zoom=alpha or 1,zoom or 1
	local N=NUM
	local n=num or NUM
	P=P or 100
	T=T or 100
	ROT = ROT or 0
	rot = rot or 0
	ROT = ROT * pd
	rot = rot * pd
	local p = math.pi*2/(N-0)*P*.01
	local t = math.pi*2/(n-0)*T*.01
	local rott = math.pi/180 + ROT
	local rotp = math.pi/180 + rot
	local cn=0
	local T = {}
	local T2 ={}
	local T3 ={} 
	for i=0,N-0 do
		T2[i+1]={}
		T3[i+1]={} 
		for j=0,n-0 do
			cn=cn+1
			--local r = r + (math.sin(j/n*math.pi*2+obj.time)+1)*100
			local x,y,z= Tor(R,r,p*i+rott,t*j+rotp) --z,x,y
			x,y,z = Shakev({x,y,z},cn,shakeparam)
			x,y,z = Field({x,y,z},pos, radius , strength, container)
			T[cn] = {x,y,z}
			T2[i+1][j+1] = {x,y,z}
			x,y,z  = Rot({x,y,z})
			T3[i+1][j+1] = {x,y,z}
			-- if info==1 and N*n<200 then obj.setfont("",10) obj.load("[ "..i.." / "..j.." ]") end
			obj.draw(x,y,z,alpha,zoom,-rx,-ry,-rz)
		end
	end
	--table.insert(T,T[1])
	T.v  = T2
	T.vr = T3
	return T
end


--輪郭抽出----------------------------------------------------------------------

--[[
	rikky_module.borderingを使って点描画をします。

	local skip 	= 10
	local z_num 	= 2
	local z_depth 	= 100
	local both	= 1
	local scale	= 1
	local size	= 10
	local col	= 0xffffff
	local alp	= 0.5
	local fig	= "円"
	local t = Border(skip, z_num, z_depth, both, size, col, alp, fig)
	
	--戻り値をFacets.YZに使うと側面に面を生成できます
	Facets.YZ(
		t.v, -- 面用のデータは戻り値に .v をつけます
		col,
		alp
	) 
]]

Border = function(
	progress,	-- 進行度(0～1)
	skip_num,	-- 輪郭の間引き
	z_num,		-- Z方向の数
	z_depth,	-- Z方向の間隔
	both,		-- 1ならZ方向を前後に
	scale,		-- 全体のスケール (1で等倍)
	figure,		-- {size,col,alp,figuretype} {1以上で図形を描画,図形色,図形透明度,図形形状,(文字列で指定)}
	spfield,	-- Spherical Field Option = {layer, radius , strength, container}
	randfx		-- RandomEffector = {layer, radius, strength,smooth}
	)

    if not tostring(progress):find("table:") then
	progress = {progress}
    end
    local gv=obj.getvalue	
    local skip,shift,shift_z,noloop = 10,0,0,1
    if not tostring(skip_num):find("table:") then
	skip	= skip_num
	offset	= 0
	shift_z = 0
	loop	= 1
    else
	skip	= skip_num[1] or 10
	offset	= skip_num[2] or 0
	shift_z = skip_num[3] or 0
	loop	= skip_num[4] or 1
    end
    progress = progress or {1}
    skip = math.max(1,skip)
    offset,ot = math.modf(offset)
    shift_z = math.floor(shift_z)
    scale = scale or 1
    z_depth = z_depth or 100
    obj.effect("領域拡張","上",2,"下",2,"左",2,"右",2)
    local w,h = obj.getpixel()
    local piw,pih = math.pi/w,math.pi/h
    local b,n,m = rikky_module.bordering()
    local nz_num = (both==1 and (-z_num) ) or 0

    local Draw = function() return 0 end
    local size, col, alp, fig = 0,0xffffff,1,"円"
    
    if (figure) and figure[1]>0 then
     	if fig=="背景" then fig="円" end
     	size, col, alp, fig = figure[1] , figure[2] or 0xffffff, figure[3] or 1, figure[4] or "円"
        local tmp = Param()
        obj.load("figure",fig ,col , size)
        Param(tmp)
	obj.setoption("billboard",3)
        Draw = function(x,y,z,Z,a,...)
            x,y,z = Rot(Scale({x,y,z},obj.zoom))
            obj.draw(x,y,z,Z,a,...) 
        end
    end

    local Field = function(v) return v[1],v[2],v[3] end
    local pos, radius , strength, container = 0,0,0,0
    if spfield then
	Field = SPField
	pos, radius , strength, container = unpack(spfield)
    end

	local r_layer,r_radius,r_strength
	if randfx then
	 	r_layer  = randfx[1]
	 	r_radius = randfx[2] or 400
		r_strength = randfx[3] or 100
		r_strength =r_strength*.01
		smooth = randfx[4] or 1
	end
    local t ={}
    local t2={}
    local t3={}
    for i = 1,#b do
	t[i]={}
        t2[i]={}
	t3[i]={}
        local p=0
	local pl=0
	local prog = math.min(1, math.max(.001, progress[i%(#progress)+1] ) )
	local num = math.ceil(n[i] * prog )
	--輪郭の頂点が3以上あった場合のみ(skipされる分もカウント)
	--面が引けずに点だけ生成されるのがあまりきれいでなかったので。
	if num>(1+skip) then
        	for j = 1,num,skip do
		local J  = (j+offset  )%n[i]+1
		--local J2 = (j+offset+skip)%n[i]+1
           	p=p+1
           	t2[i][p]={}
		t3[i][p]={}
            	local z=0
            	for k = nz_num, z_num do
               		z  = z+1
			pl = pl+1
			J  =  (J + k*shift_z )%n[i]+1
			--J2  = (J2+ k*shift_z )%n[i]+1
                	local px,py,pz = b[i][J *2-1]*scale , b[i][J *2]*scale , k * z_depth
			--local dx,dy    = b[i][J2*2-1]*scale , b[i][J2*2]*scale
			--px = px*ot + dx*(1-ot)
			--py = py*ot + dy*(1-ot)

			if randfx then
				for L =1,#r_layer do
					local l = "layer".. (r_layer[L]) .."."
					if gv(l.."x") then
						local len = Lensq( {px+obj.x+obj.ox,py+obj.y+obj.oy,pz+obj.z+obj.oz} , {gv(l.."x"), gv(l.."y"), gv(l.."z")} ) ^.5
						local st = Linear(len,0,r_radius,1,0)^smooth * r_strength
						px = px + rand(-100,100,i+p+L,p+1000)*st
						py = py + rand(-100,100,j+p+L+10,p+3000)*st
						pz = pz + rand(-100,100,k+p+L+20,p+6000)*st
					end
				end
			end

                	px, py, pz = Field({px+obj.ox, py+obj.oy, pz+obj.oz},pos, radius , strength, container)
			t3[i][p][z]= {Rot({px, py, pz})}
			px, py, pz = px-obj.ox, py-obj.oy, pz-obj.oz
                	t2[i][p][z]={px, py, pz}
                	t[i][pl]={px, py, pz}
                	--obj.putpixel(px+X,py+Y,0xffffff,1)
                	Draw(px,py,pz,.5,alp)
            	end
        	end
	end
    if (loop==1) then
       table.insert(t2[i],t2[i][1])
       table.insert(t3[i],t3[i][1])
    end
    end	
	obj.setoption("billboard",0)
    t.v  = t2
    t.vr = t3
    t.b = b
    return t
   
end


--正多面体----------------------------------------------------------------------
local G = (1+math.sqrt(5))*.5
local Gp = G*G
Regular_p = {
	--Tetrahedron
	{-1,-1,-1,  -1,1,1,  1,-1,1,  1,1,-1},

	-- Cube
	{-1,-1,-1,  -1,-1,1,  -1,1,-1,  -1,1,1,  1,-1,-1,  1, -1,1,  1,1,-1,  1 ,1,1},

	-- Octahedron
	{-1,0,0  ,0,-1,0,  0,0,-1,  1,0,0,  0,1,0,  0,0,1},

	-- Dodecahedron
	{0,-1,-Gp,  0,1,-Gp,  0,-1,Gp,  0,1,Gp,  -1,-Gp,0,  1,-Gp,0,  -1,Gp,0,  1,Gp,0,  -Gp,0,-1,  -Gp,0,1,  Gp,0,-1,  Gp,0,1,  -G,-G,-G,  -G,-G,G,  -G,G,-G, -G,G,G,  G,-G,-G,  G,-G,G,  G,G,-G,  G,G,G},

	-- Icosahedron
	{1,G,0, 1,-G,0, -1,G,0, -1,-G,0,  0,1,G,  0,1,-G,  0,-1,G,  0,-1,-G,  G,0,1,  G,0,-1,  -G,0,1,  -G,0,-1}
}
--[[
	--こんな感じ
	local v = Regular_p[3]
	local scl = 300
	local temp = {}
	for i=1,#v,3 do
		local x,y,z = v[i]*scl, v[i+1]*scl, v[i+2]*scl
		temp[math.ceil(i/3)] = {x,y,z}
		obj.draw(x,y,z)
	end
]]


--<< グリッド配置に面  >>---------------------------------------------------------------------

--[[
	各配置の戻り値 t.dを使用して面を作ります。
	というかarray[i][j][k] = {x,y,z} 形式ならなんでも

	--使用例(Dots + Facets.XY, Facets.XZ, Facets.YZ)
	local num_xyz = {10,10,3}
	lcoal div_xyz = {100,100,400}
	local t = Dots(num_xyz, div_xyz)

	local vertex = t.v 		-- Dotsの戻り値のテーブルの["v"]に配置データ
	local polycolor = 0xffffff	-- 色
	local alpha = 1			-- 透明度
	local tri = 0 			-- 1だと三角形に分割
	Facets.XY(vertex, polycolor, alpha, tri ) --XY正面
	Facets.XZ(vertex, polycolor, alpha, tri ) --XZ上面
	Facets.YZ(vertex, polycolor, alpha, tri ) --YZ側面
	
	--使用例(Cir と Facets.C)
	local rad,num = 300,12
	local vertex = Cir(rad,num)
	 -- Cirの戻り値,色,透明度
	Facets.C(vertex,0xfffff,1, .3)
	
	--使用例(Sphere + Facets.SP)
	local rad,num = 300,{12,12}
	local vertex = Sphere(rad,num)
	 --Sphereの戻り値,色,透明度
	Facets.SP(vertex,0xfffff,.5)
]]


	-- 重心座標、エッジの最大長、面法線 ---------------------------------------------------------------
	local Face = function(v)
		local a={v[1] ,v[2] ,v[3] }
		local b={v[4] ,v[5] ,v[6] }
		local c={v[7] ,v[8] ,v[9] }
		local d={v[10],v[11],v[12]}
		--[[
		local a_b = Sub(b,a)
		local b_c = Sub(c,b)
		local c_d = Sub(d,c)
		local d_a = Sub(a,d)
		local face0 = Cross(a_b,b_c)
		local face1 = Cross(c_d,d_a)
		local face = Scale( Add(face1, face0) ,.5)
		local normal = Normalize(face)
		local maxlen = math.max(Lensq(a,b),Lensq(b,c),Lensq(c,d),Lensq(d,a))
		maxlen = math.sqrt(maxlen)
		]]
		local center = {(a[1]+b[1]+c[1]+d[1])/4,(a[2]+b[2]+c[2]+d[2])/4,(a[3]+b[3]+c[3]+d[3])/4}
  		return center --,maxlen,normal
	end
--[[
local Colordiff = function(material_col,light_col,intensity)
	intensity = intensity or 1
	local mr,mg,mb = RGB(material_col)
	local lr,lg,lb = RGB(light_col)
	local m = 1/255 * intensity
	return {mr*lr*m , mg*lg*m, mb*lb*m , m}
end
]]


local Clamp = function(v,min,max)
	min,max=min or 0,max or 255
	return {
		math.max(min,math.min(max,v[1])),
		math.max(min,math.min(max,v[2])),
		math.max(min,math.min(max,v[3]))
	}
end

local Colvec = function(col)
	local r,g,b = RGB(col)
	return {r=r,g=g,b=b}
end

Blend = {}
	Blend.Multiply = function(a,b)
		-- A*B
		local ret = {a[1]*b[1],a[2]*b[2],a[3]*b[3]}
		return Clamp(ret)
	end

	Blend.ColorBurn = function(a,b)
		-- 1-(1-B)/A
		local n = {255,255,255}
		local ret = Div( Sub(n,b), a)
		ret = Sub(n,ret)
		return Clamp(ret)
	end
	
	Blend.LinearBurn = function(a,b)
		--A+B-1
		return {
			a[1]+b[1]-255,
			a[2]+b[2]-255,
			a[3]+b[3]-255
		}
	end
	
	Blend.Screen = function(a,b)
		-- 1-(1-a)*(1-b)
		local n = {255,255,255}
		local a0 = Sub(n,a)
		local b0 = Sub(n,b)
		local ret =  Sub(n,Vector.Mul(a0,b0) )
		return Clamp(ret)
	end

	Blend.ColorDodge = function(a,b)
		-- B/(1-A)
		b = Add(b,{1,1,1})
		a = Add(a,{1,1,1})
		local a0 = Sub({255,255,255}, a)
		local ret = Scale(Div(b,a0),1/255)
		return Clamp(ret)
	end
	
	Blend.LinearDodge = function(a,b,power)
		-- A+B
		b = Scale(b,power or 1)
		local ret = Add(a,b)
		return Clamp(ret)
	end

	Blend.Add = Blend.LinearDodge
	
	Blend.Lighten = function(a,b)
		-- max(A,B)
		local r = math.max(a[1],b[1])
		local g = math.max(a[2],b[2])
		local b = math.max(a[3],b[3])
		return {r,g,b}
	end
	
	Blend.Darken = function(a,b)
		-- min(A,B)
		local r = math.min(a[1],b[1])
		local g = math.min(a[2],b[2])
		local b = math.min(a[3],b[3])
		return {r,g,b}
	end
	
	Blend.Overlay = function(a,b)
		-- a * (a+ 2*b * (1-a))
	end

	--簡易ライティング-------------------------------------------------------------

Makelight = function(intensity,color,distance)
	if not LightLayer then LightLayer={} end
	intensity	= intensity or 1
	color		= color or obj.getpixel(0,0)
	distance 	= distance or 1200
	distance	= math.max(1,distance)
	local r,g,b 	= RGB(color)

	LightLayer[obj.layer] = { 
		intensity = intensity,
		pos ={obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz},
		distance = distance,
		r = r,
		g = g,
		b = b
		}
end


local Refx = function (vertex, intensity, material_col, ambient_col, LightLayer , blend, min)

	if not LightLayer then return end
	local temp={}
	local cn=0
	for k,v in pairs(LightLayer) do
		if obj.getvalue("layer"..k..".x") then
			cn=cn+1
			temp[cn]=v
		end
	LightLayer = temp
		
	end
	if #LightLayer<1 then return end

	--local light = {}
	local g=obj.getvalue

	min = min or 0.85

	local intensity_diffuse , intensity_specular ,intensity_ambient
	if not tostring(intensity):find("table:") then
		intensity_diffuse  = intensity
		intensity_specular = intensity
		intensity_ambient = 0.1
	else
		intensity_diffuse  = intensity[1]
		intensity_specular = intensity[2] or intensity_diffuse 
		intensity_ambient  = intensity[3] or 0.1
	end
		
	--[[
	if not tostring(lightlayer):find("table:") then
		local l="layer" .. (lightlayer) .."."
		if not g(l.."x") then 
			return
		else
			light[1] = {g(l.."x"),g(l.."y"),g(l.."z")}
		end
	elseif (#lightlayer>0) then
		local cn=0
		for i=1,#lightlayer do
			local l="layer" .. (lightlayer[i]) .."."
			if g(l.."x") then 
				cn = cn+1
				light[cn] = {g(l.."x"),g(l.."y"),g(l.."z")}
			end
		end
		if cn==0 then return end
	end

	if not tostring(light_col):find("table:") then
		light_col = {light_col} 
	elseif #light_col<#lightlayer then
		local num = #lightlayer - #light_col
		for i=0,num do
			table.insert(light_col,light_col[1])
		end
	end


	
	local rot = {g("rx"),g("ry"),g("rz")}
	local temp = {}
	for i=1,12,3 do
		local x,y,z = Rot({vertex[i],vertex[i+1],vertex[i+2]},rot)
		temp[i]   = x
		temp[i+1] = y
		temp[i+2] = z
	end
	vertex = temp
	temp = nil
	]]

	material_col = material_col or obj.getpixel(0,0)
	local R,G,B    = RGB(material_col)
	obj.pixeloption("type","rgb")
	obj.putpixel(0,0,R,G,B,255)
	R,G,B = R/255,G/255,B/255
	--local ambient_col = light_col[2]
	local aR,aG,aB = RGB(0xffffff)
	local aR,aG,aB = aR*R*intensity_ambient, aG*G*intensity_ambient, aB*B*intensity_ambient -- ambient

	local cam = rikky_module.camerainfo(vertex)		 -- vertexはdrawpoly用の頂点4セット
	local surface = {cam.mx, cam.my, cam.mz}		 -- 面法線
	local eye = {cam.vx, cam.vy, cam.vz}			 -- 面へのカメラ視線ベクトル
	local center = Face(vertex)				 -- 面の重心座標
	center = Add(center,{obj.x, obj.y, obj.z})		 -- 基準座標のみ追加
	local vc = Dot(eye ,surface)				 -- 面がカメラ方向に向いているか (カメラ==ライトの場合はこれだけ)
	local dR,dG,dB = 0,0,0
	local rR,rG,rB = 0,0,0
	local mR,mG,mB = R,G,B
	--obj.putpixel(0,0,R,G,B,255)
	local rt = {}
	for k,light in ipairs(LightLayer) do
		local lR,lG,lB = light.r, light.g, light.b
		lR,lG,lB = lR/255,lG/255,lB/255
		local lcol = {lR,lG,lB}
		local light_vec,dist = Normalize(Sub(light.pos,center))	 -- 面→ライトのベクトル、面との距離
		local vd0 = Dot(surface,light_vec)	 		 -- vd,面がライトに向いているか(単純なDiffuse用)
		local inter_eye = (Dot(eye,surface) * vd0)>0 and 0 or 1	 -- 光源とカメラがポリゴンで遮られているか簡易判定
		vd0 = math.max(0,math.min(1,vd0)) --* inter_eye
		local power = 1/math.max(1,(dist/light.distance)^2)		 -- 距離減衰
		local vd = vd0 * intensity_diffuse * power * light.intensity
		local vr = Refrect_N(surface ,light_vec)		 -- vr,面に反射した後のライトのベクトル
		local v2 = Dot(eye ,vr)				 	 -- v2,反射光がカメラに向いているか(反射光用)
		v2 = math.max(0,v2) * inter_eye
		v2 = Linear(v2,min,1,0,1) * intensity_specular * power
		V2_RET = v2
	
		dR,dG,dB = dR*1 + aR*0 + 255*mR*lR*vd,  dG*1 + aG*0 + 255*mG*lG*vd,  dB*1 + aB*0 + 255*mB*lB*vd
		--dR,dG,dB = math.min(255,dR),math.min(255,dG),math.min(255,dB)
		--obj.putpixel(0,0,math.min(255,dR+aR),math.min(255,dG+aG),math.min(255,dB+aB),vd*255)

		rR,rG,rB = (dR*1+rR + 255*lR)*v2, (dG*1+rG + 255*lG)*v2, (dB*1+rB + 255*lB)*v2
		rR,rG,rB = math.min(rR,255), math.min(rG,255), math.min(rB,255)
		rt[k] = {r=rR, g=rG ,b=rB, a=v2*255}				 -- スペキュラはライト個別で加算

	end

	obj.pixeloption("blend",0)
	obj.putpixel(0,0,math.min(255,dR+aR),math.min(255,dG+aG),math.min(255,dB+aB),255)

	for k,v in ipairs(rt) do
		obj.pixeloption("blend",blend or 0)
		obj.putpixel(0,0, v.r, v.g, v.b, math.min(255,v.a))   -- スペキュラ
	end

	obj.pixeloption("blend")
	obj.pixeloption("type","col")
end

--本体 ----------------------------------------------------------------------------------

Facets = {
	--正面
	XY = function(t,col,alp,tri,reflect)
		local x,y,z = 0,0,0 --Rotc() --obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz
		obj.setoption("billboard",0)
		alp = alp or 1
		tri = tri or 0
		local Ref = function() return end
		local Put = function() return end
		reflect = reflect or {}
		if #reflect>0 then Ref = Refx end
		obj.setoption("billboard",0)
		alp = alp or 1
		alpha_pat = alpha_pat or 0

		local t_col,t_alp = obj.getpixel(0,0)
		local multi_col=0
		if not col and t_alp==0 then obj.putpixel(0,0,0xffffff,1) end 
		if col then 
			if not tostring(col):find("table:") then
				obj.putpixel(0,0,col,1)
				col = {col}
			else
				Put = obj.putpixel
				
			end
		end

		alp = alp or 1
		tri = tri or 0
		--if col then obj.putpixel(0,0,col,1) end
		local xn=#t
		local yn=#t[1]
		local zn=#t[1][1]
		if (tri==0) then
			for i=1,xn-1 do
			for j=1,#t[1]-1 do
			for k=1,#t[i][j] do
				local x0,y0,z0 = t[i  ][j  ][k][1]	,t[i  ][j  ][k][2]	,t[i  ][  j][k][3]
				local x1,y1,z1 = t[i+1][j  ][k][1]	,t[i+1][j  ][k][2]	,t[i+1][j  ][k][3]
				local x2,y2,z2 = t[i+1][j+1][k][1]	,t[i+1][j+1][k][2]	,t[i+1][j+1][k][3]
				local x3,y3,z3 = t[i  ][j+1][k][1]	,t[i  ][j+1][k][2]	,t[i  ][j+1][k][3]
				Put(0,0,col[i%(#col)+1],1)
				Ref({x0+x,y0+y,z0+z, x1+x,y1+y,z1+z, x2+x,y2+y,z2+z, x3+x,y3+y,z3+z},unpack(reflect))
				obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
			end
			end
			end
		elseif (tri==1) then
			for i=1,xn-1 do
			for j=1,#t[1]-1 do
			for k=1,#t[i][j] do
				local x0,y0,z0 = t[i  ][j  ][k][1]	,t[i  ][j  ][k][2]	,t[i  ][  j][k][3]
				local x1,y1,z1 = t[i+1][j  ][k][1]	,t[i+1][j  ][k][2]	,t[i+1][j  ][k][3]
				local x2,y2,z2 = t[i+1][j+1][k][1]	,t[i+1][j+1][k][2]	,t[i+1][j+1][k][3]
				local x3,y3,z3 = t[i  ][j+1][k][1]	,t[i  ][j+1][k][2]	,t[i  ][j+1][k][3]
				if (j%2==0) then
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2, 0,0,0,0,0,0,0,0,alp)
					obj.drawpoly(x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				else
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
					obj.drawpoly(x1,y1,z1, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				end
			end
			end
			end
		end

	end,

	--上面
	XZ = function(t,col,alp,tri,reflect)
		local x,y,z = 0,0,0 -- Rotc() --obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz
		obj.setoption("billboard",0)
		alp = alp or 1
		tri = tri or 0
		local Ref = function() return end
		local Put = function() return end
		reflect = reflect or {}
		if #reflect>0 then Ref = Refx end
		obj.setoption("billboard",0)
		alp = alp or 1
		alpha_pat = alpha_pat or 0

		local t_col,t_alp=obj.getpixel(0,0)
		local multi_col=0
		if not col and t_alp==0 then obj.putpixel(0,0,0xffffff,1) end 
		if col then 
			if not tostring(col):find("table:") then
				obj.putpixel(0,0,col,1)
				col = {col}
			else
				Put = obj.putpixel
				
			end
		end

		alp = alp or 1
		tri = tri or 0
		--if col then obj.putpixel(0,0,col,1) end
		local xn=#t
		local yn=#t[1]
		local zn=#t[1][1]
		if (tri==0) then
			for i=1 ,xn-1	do
			for j=1 ,yn	do
			for k=zn,2,-1	do --z軸は奥から引く
				local x0,y0,z0 = t[i  ][j][k  ][1]	,t[i  ][j][k  ][2]	,t[i  ][j][k  ][3]
				local x1,y1,z1 = t[i+1][j][k  ][1]	,t[i+1][j][k  ][2]	,t[i+1][j][k  ][3]
				local x2,y2,z2 = t[i+1][j][k-1][1]	,t[i+1][j][k-1][2]	,t[i+1][j][k-1][3]
				local x3,y3,z3 = t[i  ][j][k-1][1]	,t[i  ][j][k-1][2]	,t[i  ][j][k-1][3]
				Put(0,0,col[i%(#col)+1],1)
				Ref({x0+x,y0+y,z0+z, x1+x,y1+y,z1+z, x2+x,y2+y,z2+z, x3+x,y3+y,z3+z},unpack(reflect))
				obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
			end
			end
			end
		elseif (tri==1) then
			for i=1	,xn-1	do
			for j=1	,yn	do
			for k=zn,2,-1	do --z軸は奥から引く
				local x0,y0,z0 = t[i  ][j][k  ][1]	,t[i  ][j][k  ][2]	,t[i  ][j][k  ][3]
				local x1,y1,z1 = t[i+1][j][k  ][1]	,t[i+1][j][k  ][2]	,t[i+1][j][k  ][3]
				local x2,y2,z2 = t[i+1][j][k-1][1]	,t[i+1][j][k-1][2]	,t[i+1][j][k-1][3]
				local x3,y3,z3 = t[i  ][j][k-1][1]	,t[i  ][j][k-1][2]	,t[i  ][j][k-1][3]
				--local cam = rikky_module.camerainfo({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3})
				if (k%2==0) then
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2, 0,0,0,0,0,0,0,0,alp)
					obj.drawpoly(x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				else
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
					obj.drawpoly(x1,y1,z1, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				end
			end
			end
			end
		end
	end,

	--側面
	YZ = function(t,col,alp,tri,reflect)
		local x,y,z = 0,0,0 --Rotc() --obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz
		obj.setoption("billboard",0)
		alp = alp or 1
		tri = tri or 0
		local Ref = function() return end
		local Put = function() return end
		reflect = reflect or {}
		if #reflect>0 then Ref = Refx end
		obj.setoption("billboard",0)
		alp = alp or 1
		alpha_pat = alpha_pat or 0

		local t_col,t_alp=obj.getpixel(0,0)
		local multi_col=0
		if not col and t_alp==0 then obj.putpixel(0,0,0xffffff,1) end 
		if col then 
			if not tostring(col):find("table:") then
				obj.putpixel(0,0,col,1)
				col = {col}
			else
				Put = obj.putpixel
				
			end
		end
		local xn=#t
		--local yn=#t[1]
		--local zn=#t[1][1]
		if (tri==0) then
			for i=1,xn		do
			for j=1,#t[i]-1		do
			for k=#t[i][j],2,-1	do --z軸は奥から引く
				local x0,y0,z0 = t[i][j  ][k  ][1]	,t[i][j  ][k  ][2]	,t[i][j  ][k  ][3]
				local x1,y1,z1 = t[i][j  ][k-1][1]	,t[i][j  ][k-1][2]	,t[i][j  ][k-1][3]
				local x2,y2,z2 = t[i][j+1][k-1][1]	,t[i][j+1][k-1][2]	,t[i][j+1][k-1][3]
				local x3,y3,z3 = t[i][j+1][k  ][1]	,t[i][j+1][k  ][2]	,t[i][j+1][k  ][3]
				--local cam = rikky_module.camerainfo({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3})
				Put(0,0,col[i%(#col)+1],1)
				
					Ref({x0,y0+y,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3},unpack(reflect))
				obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
			end
			end
			end
		elseif (tri==1) then
			for i=1,xn		do
			for j=1,#t[i]-1		do
			for k=#t[i][j],2,-1	do --z軸は奥から引く
				local x0,y0,z0 = t[i][j  ][k  ][1]	,t[i][j  ][k  ][2]	,t[i][j  ][k  ][3]
				local x1,y1,z1 = t[i][j  ][k-1][1]	,t[i][j  ][k-1][2]	,t[i][j  ][k-1][3]
				local x2,y2,z2 = t[i][j+1][k-1][1]	,t[i][j+1][k-1][2]	,t[i][j+1][k-1][3]
				local x3,y3,z3 = t[i][j+1][k  ][1]	,t[i][j+1][k  ][2]	,t[i][j+1][k  ][3]
				Put(0,0,col[i%(#col)+1],1)
				if (i%2==0) then
						Ref({x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2},unpack(reflect))
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2, 0,0,0,0,0,0,0,0,alp)
						Ref({x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3},unpack(reflect))
					obj.drawpoly(x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				else
						Ref({x1,y1,z1, x3,y3,z3, x0,y0,z0, x0,y0,z0},unpack(reflect))
					obj.drawpoly(x1,y1,z1, x3,y3,z3, x0,y0,z0, x0,y0,z0, 0,0,0,0,0,0,0,0,alp)
						Ref({x1,y1,z1, x1,y1,z1, x2,y2,z2, x3,y3,z3},unpack(reflect))
					obj.drawpoly(x1,y1,z1, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alp)
				end
			end
			end
			end
		end
	end,

	-- Cirで作った円形配置用
	C = function(t,col,alp,innner_scale,alpha_pat)
		obj.setoption("billboard",0)
		alp = alp or 1
		local sc = innner_scale or 1
		sc = 1-sc
		alpha_pat = alpha_pat or 0
		if col then obj.putpixel(0,0,col,1) end
		for i=1,#t-1 do
			local alpha=alp
			if alpha_pat==1 then
				alpha = i%2
			elseif alpha_pat==2 then
				alpha = rand(10,100,i,i)*.01
			elseif alpha_pat==3 then
				alpha = (1+math.sin(i/3.1415))*.5
			end
			
			local x0,y0,z0 = t[i  ][1],	t[i  ][2],	0
			local x1,y1,z1 = t[i+1][1],	t[i+1][2],	0
			local x2,y2,z2 = x1*sc,		y1*sc,		0
			local x3,y3,z3 = x0*sc,		y0*sc,		0
			--local cam = rikky_module.camerainfo({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3})
			obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha*alp)
		end
	end,

	-- Sphereで作った球状配置用
	SP = function(t,col,alp,tri,alpha_pat,reflect,fog)
		local Ref = function() return end
		local Put = function() return end
		reflect = reflect or {}
		if #reflect>0 then Ref = Refx end
		obj.setoption("billboard",0)
		alp = alp or 1
		tri = tri or 0
		alpha_pat = alpha_pat or 0

		local t_col,t_alp=obj.getpixel(0,0)
		local multi_col=0
		if not col and t_alp==0 then obj.putpixel(0,0,0xffffff,1) end 
		if col then 
			if not tostring(col):find("table:") then
				obj.putpixel(0,0,col,1)
				col = {col}
			else
				Put = obj.putpixel
				
			end
		end

		--if col then obj.putpixel(0,0,col,1) end
		table.insert(t,t[1])
		local wn,hn = #t, #t[1]
		local cn=0
		if (tri==0) then
			for i=1,wn-1 do
				for j=1,hn-1 do
					cn=cn+1
					local alpha=alp
					if alpha_pat==1 then
						alpha = i%2
					elseif alpha_pat==2 then
						alpha = j%2
					elseif alpha_pat==3 then
						alpha = rand(1,0,i*j,i)
					elseif alpha_pat==4 then
						alpha = rand(10,100,i*j,i)*.01
					elseif alpha_pat>=5 then
						local rep = (alpha_pat-2)
						local spd  =obj.time * (alpha_pat-4)
						alpha = (1+math.sin(i*rep + spd))*.5 * (1+math.sin(j*rep + spd))*.5
					end
					local x0,y0,z0 = t[i  ][j  ][1],	t[i  ][j  ][2],	 t[i  ][j  ][3]
					local x1,y1,z1 = t[i+1][j  ][1],	t[i+1][j  ][2],	 t[i+1][j  ][3]
					local x2,y2,z2 = t[i+1][j+1][1],	t[i+1][j+1][2],	 t[i+1][j+1][3]	
					local x3,y3,z3 = t[i  ][j+1][1],	t[i  ][j+1][2],	 t[i  ][j+1][3]
					if t[i][j].pix then
						Put(0,0,t[i][j].pix.col,1)
					else
						Put(0,0,col[cn%(#col)+1],1)
					end
					Put(0,0,col[cn%(#col)+1],1)
					Ref({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3},unpack(reflect))
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha*alp)
				end
			end
		else
			for i=1,wn-1 do
				for j=1,hn-1 do
					cn=cn+1
					local alpha=alp
					if alpha_pat==1 then
						alpha = i%2
					elseif alpha_pat==2 then
						alpha = j%2
					elseif alpha_pat==3 then
						alpha = rand(1,0,i*j,i)
					elseif alpha_pat==4 then
						alpha = rand(10,100,i*j,i)*.01
					elseif alpha_pat>=5 then
						local rep = (alpha_pat-2)
						local spd  =obj.time * (alpha_pat-4)
						alpha = (1+math.sin(i*rep + spd))*.5 * (1+math.sin(j*rep + spd))*.5
					end
					local x0,y0,z0 = t[i  ][j  ][1],	t[i  ][j  ][2],	 t[i  ][j  ][3]
					local x1,y1,z1 = t[i+1][j  ][1],	t[i+1][j  ][2],	 t[i+1][j  ][3]
					local x2,y2,z2 = t[i+1][j+1][1],	t[i+1][j+1][2],	 t[i+1][j+1][3]	
					local x3,y3,z3 = t[i  ][j+1][1],	t[i  ][j+1][2],	 t[i  ][j+1][3]
					if t[i][j].pix then
						Put(0,0,t[i][j].pix.col,1)
					else
						Put(0,0,col[cn%(#col)+1],1)
					end

					Ref({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3},unpack(reflect))
					obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2, 0,0,0,0,0,0,0,0,alpha*alp)
					obj.drawpoly(x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha*alp)
				end
			end
		end
	end,

	-- Torus配置の戻り値用 (t = t.v)
	T = function(t,col,alp,alpha_pat,reflect)
		local x,y,z = 0,0,0 --Rotc()
		obj.setoption("billboard",0)
		alp = alp or 1
		tri = tri or 0
		local Ref = function() return end
		local Put = function() return end
		reflect = reflect or {}
		if #reflect>0 then Ref = Refx end
		obj.setoption("billboard",0)
		alp = alp or 1
		alpha_pat = alpha_pat or 0

		local t_col,t_alp=obj.getpixel(0,0)
		local multi_col=0
		if not col and t_alp==0 then obj.putpixel(0,0,0xffffff,1) end 
		if col then 
			if not tostring(col):find("table:") then
				obj.putpixel(0,0,col,1)
				col = {col}
			else
				Put = obj.putpixel
				
			end
		end

		alp = alp or 1
		tri = tri or 0
		alpha_pat = alpha_pat or 0
		local wn,hn = #t, #t[1]
		for i=1,wn-1 do
			for j=1,hn-1 do
				local alpha=alp
				if alpha_pat==1 then
					alpha = i%2
				elseif alpha_pat==2 then
					alpha = j%2
				elseif alpha_pat==3 then
					alpha = rand(1,0,i*j,i)
				elseif alpha_pat==4 then
					alpha = rand(10,100,i*j,i)*.01
				elseif alpha_pat>=5 then
					local rep = (alpha_pat-2)
					local spd  =obj.time * (alpha_pat-4)
					alpha = (1+math.sin(i*rep + spd))*.5 * (1+math.sin(j*rep + spd))*.5
				end
				
				local x0,y0,z0 = t[i  ][j  ][1], t[i  ][j  ][2], t[i  ][j  ][3]
				local x1,y1,z1 = t[i+1][j  ][1], t[i+1][j  ][2], t[i+1][j  ][3]
				local x2,y2,z2 = t[i+1][j+1][1], t[i+1][j+1][2], t[i+1][j+1][3]
				local x3,y3,z3 = t[i  ][j+1][1], t[i  ][j+1][2], t[i  ][j+1][3]
				Put(0,0,col[i%(#col)+1],1)
				Ref({x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3},unpack(reflect))
				obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0, 0,0, 0,0, 0,0, alpha*alp)
			end
		end
	end
}


Facets.XYZ = function(t,col,alp,tri,ref)
	Facets.XY(t,col,alp,tri,ref)
	Facets.YZ(t,col,alp,tri,ref)
	Facets.XZ(t,col,alp,tri,ref)
end



--三角形生成(底辺wと高さhで指定)----------------------------------------------------------------------
Tri = function(w,h,col,x,y,a,rot,line_w,line_col,line_alp)
	if h==nil then h=math.sqrt(w*w*math.pi)/2 end
	x,y,a,rot=x or 0,y or 0,a or 1, rot or 0
	if col then obj.putpixel(0,0,col,1) end
	line_w,line_col,line_alp=line_w or 0,line_col or col or 0xffffff,line_alp or 1
	local R = function(x,y,r)
		r=r*math.pi/180
		return x*math.cos(r)-y*math.sin(r), x*math.sin(r)+y*math.cos(r)
	end
	w=w/2
	h=h/2
	local x0,y0=R(0,-h,rot)
	local x1,y1=R(w,h,rot)
	local x2,y2=R(-w,h,rot)
	local cx,cy = (x0+x1+x2)/3,(y0+y1+y2)/3
	x = x-cx
	y = y-cy
	x0,y0=x0+x,y0+y
	x1,y1=x1+x,y1+y
	x2,y2=x2+x,y2+y
	obj.drawpoly(x0,y0,0, x1,y1,0, x2,y2,0, x2,y2,0, 0,0,0,0,0,0,0,0,a or 1)
	if line_w>0 then
 		P_line({x0,y0},{x1,y1},line_w,line_col,line_alp)
 		P_line({x1,y1},{x2,y2},line_w,line_col,line_alp)
 		P_line({x2,y2},{x0,y0},line_w,line_col,line_alp)
	end
end



-- rikky_module.image、仮想バッファのコピーを追加（ただし仮想バッファ自体はつぶされる） ----------------------------

R_img = function(
	mode,	-- 文字で指定 "g","r","w","m"
	id,	-- 保存または読込先のid
	id1,	-- modeが "m" の場合の合成画像id
	x,	-- modeが "m" の場合の合成画像の位置x
	y	-- modeが "m" の場合の合成画像の位置y
	)

--[[

	modeが"tmp"の場合、idに仮想バッファの画像をコピーします。元画像はそのまま。
	copybuffer("cache:"..id,"tmp")の代用。
	
	local id = R_img("g")	-- 空きID検索
	R_img("tmp",id)		-- idに仮想バッファを保存	
	
	仮想バッファ自体はもともとロードしていたイメージに書き換わってしまっているので
	copybufferやloadでは元の仮想バッファのイメージは呼び出せません。
	
	rikky_module.image("r",id) または R_img("r",id)を使用します。	

]]

	--迂闊に引数なしで呼ぶと落ちるので 空きidを返すように
	if not mode or mode=="g" then
		return rikky_module.image("g")
	end
	--idはマイナスを受け付けないそうなので
	if (id) then 
		id = id < 0 and 0 or id 
	end
	if (id1) then 
		id1 = id1 < 0 and 0 or id1 
	end	
	mode = string.lower(mode)
	if (mode == "tmp") then
		if not id then return end

		-- ロード中のオリジナル画像をorg_idに退避
		local org_id = rikky_module.image("g")
		rikky_module.image("w",org_id)

		-- 目的である仮想バッファをロードして任意のidへ
		obj.copybuffer("obj","tmp")
		rikky_module.image("w",id)

		-- 退避させたオリジナル画像をロードして仮想バッファへコピー。
		-- 無駄っぽいけれどcopybuffer("cache:") を使わないために。
		-- 一度仮想バッファへ退避しないと元オブジェクトを維持しつつクリアできない(drawをしないと消える)
		rikky_module.image("r+",org_id)
		obj.copybuffer("tmp","obj")

		-- 仮想バッファの呼び出し、保存したオリジナル画像をクリア
		obj.copybuffer("obj","tmp")
		rikky_module.image("c",org_id)
		
	else
		rikky_module.image(mode,id,id1,x,y)
	end
end


--rikky_module.imageでテキストを画像にキャッシュ----------------------------------------------------------------------------------------------------

Cache_text = function(t,id,clear)

--[[
	t = {"ここに文字",...}
	テーブルの要素の添え字に s,eを使うとsで文字の最初、eで最後に指定の文字を一括につける
	n を使用すると番号を振って表示(空文字で番号のみ)
	t = {" AA"," BB"}
	t.s = "["
	t.e = "]"
	t.n = "No."
	結果 保存される画像　[No.1 AA] , [No.2 BB],  [No.3]
	idを指定した場合の保存idは id+1,id+2,id+3,...	
	idが指定されないとidは左から レイヤー番号*100 + 1, レイヤー番号*100 + 2,...
	保存した画像を呼ぶ場合は
	for i=0,#t-1 do
		rikky_module.image("r",id + (i+1) )
	end
	のようにして使用
	
	戻り値はこのスクリプトで保存された画像のイメージのidをテーブル
]]

	require("rikky_module")
	id = id or (obj.layer*100)
	local s = t.s or ""
	local e = t.e or ""
	local n = ""
	local used_id = {}
	local is_table = tostring(t):find("table:")

	if (is_table) then
		num=#t
	else
		num=t
	end

	if (obj.time==obj.totaltime) or (clear==1)  then
		for i=1,num do
			rikky_module.image("c",id+i)
		end
		return "nil"
	elseif (obj.frame<3) then
		for i=1,num do
			if (t.n) then n = (i-1) .. (t.n) end
			local txt = (is_table and t[i]) or ("")
			obj.load(s .. n .. txt .. e)
			rikky_module.image("w",id+i)
			used_id[i]=id+i
		end
	end
	return used_id
end

--球形押し出し----------------------------------------------------------------------------------------------------
--戻り値 x,y,z,d　(dは影響度を正規化したもの)
SPField = function(
	pos,		-- オブジェクトの座標 {x,y,z}
	layer,		-- 押し出しフィードの座標{x,y,z} または レイヤー番号
	radius,		-- 押し出し半径
	strength,	-- 押し出し強度
	container,	-- 1だと押し込みフィールドモード
	randomize	-- フィールドの影響をばらけさせる {amount,seed,bias}
	)
	strength = strength or 100
	strength = strength*.01
	local gv=obj.getvalue
	local pos=pos or {obj.ox,obj.oy,obj.oz}
	local ox,oy,oz = unpack(pos)
	local xx,yy,zz = (ox+obj.x),(oy+obj.y),(oz+obj.z)
	local x,y,z = 0,0,0
	radius = radius or 200
	if (layer==nil) then
		return ox,oy,oz,1
	elseif tostring(layer):find("table:") then
		x,y,z = unpack(layer)
	else
		local L = "layer"..layer
		if gv(L..".x")==nil then
			return ox,oy,oz,1
		else
			x,y,z = gv(L..".x"),gv(L..".y"),gv(L..".z")
			radius = radius + gv(L..".zoom")/2
		end
  	end
	
	local xa,ya,za = xx-x, yy-y, zz-z
	local xb,yb,zb = x-xx, y-yy, z-zz

	if radius<0 then
  	 	xa,ya,za,xb,yb,zb = xb,yb,zb,xa,ya,za
	end
  	radius = math.abs(radius)
  	local l = (xa*xa + ya*ya + za*za)^.5
	local nx,ny,nz = xb/l,yb/l,zb/l
	local d = (l<radius) and (1 - l/radius)  or 0
	d = d * strength
	
	if (container==1) then
		if (randomize) then
			local amount,seed,bias = randomize[1]*.01,randomize[2] or 1,randomize[3] or 0
			local bfx = (1-math.abs(bias))*.5
			if (l>radius) then
	    			ox = (xb-nx*radius * (bias + rand(0,100,seed,100)*amount) * bfx ) * strength + ox
	    			oy = (yb-ny*radius * (bias + rand(0,100,seed,200)*amount) * bfx ) * strength + oy
	   			oz = (zb-nz*radius * (bias + rand(0,100,seed,300)*amount) * bfx ) * strength + oz
	  		end
		else
			if (l>radius) then
	    			ox = (xb-nx*radius) * strength + ox
	    			oy = (yb-ny*radius) * strength + oy
	   			oz = (zb-nz*radius) * strength + oz
	  		end
		end
		return ox,oy,oz,d
	else
		if (l<radius) then
	    		ox = (xb-nx*radius) * strength + ox
	    		oy = (yb-ny*radius) * strength + oy
	    		oz = (zb-nz*radius) * strength + oz
	  	end
		return ox,oy,oz,d
	end

end
------------------------------------------------------------------------------------------------------
-- falseを返すようにするとrequireする毎に更新される
-- return false
