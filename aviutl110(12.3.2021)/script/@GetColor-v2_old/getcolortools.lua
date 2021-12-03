
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

-- カメラ平面作成---------------------------------------------------
-- 引数なしでカメラ目標の距離に平面、offsetがあるとカメラから1024+offsetの位置
local Norm_cameraplane = function(offset)
    local c = obj.getoption("camera_param")
    local n,d = Normalize( Sub( {c.tx,c.ty,c.tz}, {c.x,c.y,c.z} ))
    if offset then d=1024+offset end
    return {n[1],n[2],n[3],d=d}
end

-- 3点a,b,cのポリゴン法線,重心 ----------------------------
local Norm_surface = function(a,b,c)
    local AB=Sub(b,a)
    local BC=Sub(c,b)
    local center = {(a[1]+b[1]+c[1])/3,(a[2]+b[2]+c[2])/3,(a[3]+b[3]+c[3])/3}
    local N = Normalize(Cross(AB,BC))
    return N,center
end

-- 点Aと平面上の最近点(A=座標 P=平面上の点 N=平面の法線 )----------------------------
-- 平面上の点 p = Scale(N,d)
local Pos_p_on_plane = function(a,p,N)
    local PA = {a[1]-p[1],a[2]-p[2], a[3]-p[3]}
    local d = (N[1]*PA[1] + N[2]*PA[2] + N[3]*PA[3])
    return {a[1]-(N[1]*d), a[2]-(N[2]*d), a[3]-(N[3]*d)}
end

-- 点pと面abc上の最近点------------------------------------
local Pos_p_on_poly = function(p,a,b,c)
  local AB = Sub(b,a)
  local BC = Sub(c,b)
  local N = Normalize(Cross(AB,BC))
  --local PA = Sub(A,a)
  --local d = Dot( N, PA )
  --return Sub(A,Scale(N,d)),N
  return Pos_p_on_plane(p,a,N)
end

-- 点Pと線abの距離------------------------------------
local Dist_p_line = function(P,a,b)
  local ab,aP = Sub(b,a), Sub(P,a)
  local cross   = Cross(ab,aP)
  local len   = Length(ab)
  if len==0 then len=eps end -- 0除算でnanが出るので
  local dst   = Length(cross) / len
  return dst
end

-- 点Pと線分abの距離------------------------------------
local Dist_p_segment = function(P,a,b)
    local AB = Sub(b,a)
  if	 ( Dot(AB,Sub(P,a) ) < eps ) then
    return Length(Sub(a,P))
  elseif ( Dot(AB,Sub(P,b) ) > eps ) then
    return Length(Sub(b,P))
  else
    return Dist_p_line(P,a,b)
  end
end


-- 点Pと線ab上の最近点------------------------------------
local Pos_p_on_line=function(P,a,b)
    local AB = Sub(b,a)
    local N = Normalize(AB) 	--線abの単位ベクトル
    local D = Dot(N, Sub(P,a))	--aPベクトルと内積
    return Add(a,Scale(N,D))
end

-- 点Pと線分ab上の最近点------------------------------------
local Pos_p_on_segment=function(P,a,b)
    local AB = Sub(b,a)
    if	( Dot(AB,Sub(P,a) ) < eps ) then
        return a
    elseif	( Dot(AB,Sub(P,b) ) > eps ) then
        return b
    else
        return Pos_p_on_line(P,a,b)
    end
end


-- 線分ABと平面の交点------------------------------------
--Intersect_plane_Line
local Pos_plane_intersection_segment = function(

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
local Pos_intersection_2line = function(A,B,C,D)

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
local Vec_reflect=function(surface_normal,incidence_vector)
    local s,v = surface_normal,incidence_vector
    local t = -(s[1]*v[1] + s[2]*v[2] + s[3]*v[3])/(s[1]*s[1] + s[2]*s[2] + s[3]*s[3])
  return {v[1]+(t*s[1]*2), v[2]+(t*s[2]*2), v[3]+(t*s[3]*2)}
end


--回転行列----------------------------------------------------------------------

local Rot = function(v,r)
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
local Rotc = function(pos,anc,rot)
	pos = pos or {obj.ox, obj.oy, obj.oz}
	anc = anc or {obj.cx, obj.cy, obj.cz}
	rot = rot or {obj.rx, obj.ry, obj.rz}
	local zoom = obj.getvalue("zoom")*.01
	local ox,oy,oz = pos[1],pos[2],pos[3]
	local cx,cy,cz = anc[1]*zoom,anc[2]*zoom,anc[3]*zoom
	local x,y,z = Rot({ox-cx, oy-cy, oz-cz},rot)
	return  x+cx, y+cy, z+cz
end

--vをロール(r)、ピッチ(p)、ヨー(Y)で回転(回転角はラジアン)
local Rot_rpy = function(v,r,p,Y)
    local Sin = {r=math.sin(r),p=math.sin(p),y=math.sin(Y)}
    local Cos = {r=math.cos(r),p=math.cos(p),y=math.cos(Y)}
    local x,y,z = v[1],v[2],v[3]
    local x0 = x*Cos.p*Cos.r + y*( Sin.y*Sin.p*Cos.r - Cos.y*Sin.r ) + z*( Sin.y*Sin.r + Cos.y*Sin.p*Cos.r )
    local y0 = x*Cos.p*Sin.r + y*( Sin.y*Sin.p*Sin.r + Cos.y*Cos.r ) + z*(-Sin.y*Cos.r + Cos.y*Sin.p*Sin.r )
    local z0 = -x*(Sin.p) + y*(Sin.y*Cos.p) + z*(Cos.y*Cos.p)
    return x0,y0,z0
end

--任意単位ベクトル周り回転(回転角はラジアン)
local Rot_v = function(v,r,scale)
    local x,y,z = v[1],v[2],v[3]
    scale = scale or 1
    local cosr = math.cos(r)
    local sinr = math.sin(r)
    local x0 = x*x*(1-cosr)+cosr   + x*y*(1-cosr)-z*sinr + z*x*(1-cosr)+y*sinr
    local y0 = x*y*(1-cosr)+z*sinr + y*y*(1-cosr)+cosr   + y*z*(1-cosr)-x*sinr
    local z0 = z*x*(1-cosr)-y*sinr + y*z*(1-cosr)+x*sinr + z*z*(1-cosr)+cosr
    return x0*scale,y0*scale,z0*scale
end

--球形配置coordinates r = radius, p = math.pi*2/hn*i + mathpi/2 , t = math.pi/hn*j
local Co_spherical = function(r,p,t)
    local x,y,z = math.sin(t)*math.cos(p), math.sin(t)*math.sin(p), math.cos(t)
    return y*r,-z*r,x*r
end

--トーラス配置,R=Radius, r=raidus, p = math.pi*2/hn*i + math.pi/2 , t = math.pi/hn*j
local Co_torus = function(R,r,p,t)
    local x = R*math.cos(t) + r*math.cos(p)*math.cos(t)
    local y = R*math.sin(t) + r*math.cos(p)*math.sin(t)
    local z = r*math.sin(p)
    return x,y,z
end
--パラボライダル u=(0~ ) ,v=(0~ ),r=(0 <= math.pi*2)
local Co_Paraboloidal=function(u,v,r)
    local x = u*v*math.cos(r)
    local y = u*v*math.sin(r)
    local z = (u*u-v*v)/2
    return x,y,z
end

-- Prolate spheroidal coordinates xi=[0,huge) ,n=[0,math.pi), p=[0,math.pi*2)
local Co_prolate_spheroidal=function(xi,n,p,a)
    a = a or 1
    local x = a * math.sinh(xi)*math.sin(n)*math.cos(p)
    local y = a * math.sinh(xi)*math.sin(n)*math.sin(p)
    local z = a * math.cosh(xi)*math.cos(n)
    return x,y,z
end

--他のレイヤー座標をまとめる---------------------------------------------------

local GL = function(...)
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

local GL2 = function(...)
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

local Param = function(param)
	if tostring(param):find("table:") then
		obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.cx,obj.cy,obj.cz,obj.aspect=unpack(param)
	elseif param==0 then
		obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.aspect,obj.cx,obj.cy,obj.cz = 0,0,0,1,1,0,0,0,0
	else
		return {obj.ox,obj.oy,obj.oz,obj.zoom,obj.alpha,obj.cx,obj.cy,obj.cz,obj.aspect}
	end
end

local Camparam = function()
	local c = obj.getoption("camera_param")
	c.pos = {c.x,c.y,c.z}
	c.tgt = {c.tx,c.ty,c.tz}
	c.eye = {c.tx-c.x,c.ty-c.y,c.tz-c.z}
	c.l = math.sqrt(c.eye[1]*c.eye[1] + c.eye[2]*c.eye[2] + c.eye[3]*c.eye[3])
	if c.l==0 then
		c.n = {0,0,-1}
	else
		c.n = {c.eye[1]/c.l, c.eye[2]/c.l, c.eye[3]/c.l}
	end
	return c
end

-- linear変換
local Linear = function(t, t1, t2, v1, v2, nolimit )
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

local TBL = function(v,n)
    --vがnilなら0を(nを指定すればn)、
    --テーブルならobj.indexで振り分け(要素がobj.numより少なければループ)
    --それ以外ならそのまま返します。
    if n==nil then n=0 end
    local V=tostring(v)
    	if v==nil then
    	return n
    	elseif not string.find(V,"table:")then
    	return v
    	else
    	return v[obj.index % (#v) + 1]
    	end
end
local Swap_table = function(t,a,b)
    if t==0 then return a end
    if t==1 then return b end
    local t1 = 1-t
    return {
        a[1]*t1+b[1]*t,
        a[2]*t1+b[2]*t,
        a[3]*t1+b[3]*t
    }
end

local Swap_col=function(t,cola,colb)
    local a={RGB(cola)}
    local b={RGB(colb)}
    local c = Swap_table(t,a,b)
    return RGB(c[1],c[2],c[3]),c
end

local Swap_i_table = function(t,...)
	t=t-1
	local i,d = math.modf(t)
	local len = select("#",...)
    if d==0 then return select(i%len+1,...) end
	local a = select(i%len+1,...)
	local b = select((i+1)%len+1,...)
    return {
        a[1]*(1-d) + b[1]*d,
        a[2]*(1-d) + b[2]*d,
        a[3]*(1-d) + b[3]*d
    }
end

local Progress = function(progress, wait, num, ease)
    --[[
    トラックバーで動作するTAの動き
    numで指定した長さのテーブルを返します。
    テーブルの中身は0~1の範囲の値で
    progressが0~1に動作するとテーブルの中身も0~1に変化します
    waitでテーブルごとにタイミングをずらせます。
    easeはeasing.luaが呼び出せる場合指定すると使えます。
    easeは番号で、または文字列でイージング名を直接指定できます。
    ]]
    wait = wait or 0.1
    num = num or 10
    local dulation = 1 + wait*(num-1)
    local t = progress * dulation
    local T = {}
    if progress==0 then
        for i=1,num do
            T[i]=0
        end
        return T
    elseif progress==1 then
        for i=1,num do
            T[i]=1
        end
        return T
    end
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
local Shuffle = function(num,order,seed,scl)

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
            index[k]=i
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

  --点p0と点p1を結ぶ直線(2D)----------------------------------------------------------------------

  local Draw_line = function(
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

  local Draw_line3D = function(
  	p0,	--座標 {x,y,z}
  	p1,	--座標 {x,y,z}
  	width,	--[線幅]
  	col,	--[色]
  	alp,	--[透明度]
  	st,	--[消滅開始距離]
  	va,	--[消滅までのフェード範囲]
  	nst,	--[消滅開始距離(st以下)]
  	nva,	--[消滅までの範囲(nst以下～0まで)]
  	t,	--[0~1で線を伸ばす]
    cam_param
  	)
  	width = width or 1
  	alp = alp or 1
  	st,va = st or 500,va or 1000
  	nst,nva = nst or 1 ,nva or 0
  	t = t or 1
  	t = math.max(0,math.min(1,t))
  	p0[3] = p0[3] or 0
  	p1[3] = p1[3] or 0
  	local a=Sub(p1,p0)
  	local len=Length(a)
  	if len>(st+va) then return 0,p1 end
  	if len<(nst-nva) then return 0,p1 end
  	local lin = Linear(len,st,va,1,0)
  	lin = lin*Linear(len,nva,nst ,0,1)
  	width = width*lin
  	alp = alp*lin
  	if col then obj.putpixel(0,0,col,1) end
  	local c = cam_param or obj.getoption("camera_param")
  	local b={c.x-p0[1], c.y-p0[2], c.z-p0[3]}
  	local n = Cross(a,b)
  	local l = Length(n)
  	local nx,ny,nz = (n[1]/l)*width*.5, (n[2]/l)*width*.5 ,(n[3]/l)*width*.5
  	p1 = Add( Scale(p1,t), Scale(p0,1-t))
  	obj.drawpoly(
  		p0[1]-nx,p0[2]-ny,p0[3]-nz,
  		p1[1]-nx,p1[2]-ny,p1[3]-nz,
  		p1[1]+nx,p1[2]+ny,p1[3]+nz,
  		p0[1]+nx,p0[2]+ny,p0[3]+nz,
  		0,0,0,0,0,0,0,0,alp
  		)
  	return 1
  end

  local Draw_line3D2 = function(
  	p0,	--座標 {x,y,z}
  	p1,	--座標 {x,y,z}
  	width,	--[線幅]
  	col,	--[色]
    col2,
  	alp,	--[透明度]
    alp2,
  	st,	--[消滅開始距離]
  	va,	--[消滅までのフェード範囲]
  	nst,	--[消滅開始距離(st以下)]
  	nva,	--[消滅までの範囲(nst以下～0まで)]
  	t,	--[0~1で線を伸ばす]
    cam_param
  	)
  	width = width or 1
  	st,va = st or 500,va or 1000
  	nst,nva = nst or 1 ,nva or 0
  	t = t or 1
  	t = math.max(0,math.min(1,t))
  	p0[3] = p0[3] or 0
  	p1[3] = p1[3] or 0
  	local a=Sub(p1,p0)
  	local len=Length(a)
  	if len>(st+va) then return 0,p1 end
  	if len<(nst-nva) then return 0,p1 end
    alp2 = alp2 or alp1
    col2 = col2 or col
  	local lin = Linear(len,st,va,1,0)
  	lin = lin*Linear(len,nva,nst ,0,1)
  	width = width*lin
  	local alpha = math.max(alp,alp2)--lin
  	obj.putpixel(0,0,col,alp)
  	obj.putpixel(1,0,col2,alp2)
    obj.putpixel(2,0,col2,alp2)
  	local c = cam_param or obj.getoption("camera_param")
  	local b={c.x-p0[1], c.y-p0[2], c.z-p0[3]}
  	local n = Cross(a,b)
  	local l = Length(n)
  	local nx,ny,nz = (n[1]/l)*width*.5, (n[2]/l)*width*.5 ,(n[3]/l)*width*.5
  	p1 = Add( Scale(p1,t), Scale(p0,1-t))
  	obj.drawpoly(
  		p0[1]-nx,p0[2]-ny,p0[3]-nz,
  		p1[1]-nx,p1[2]-ny,p1[3]-nz,
  		p1[1]+nx,p1[2]+ny,p1[3]+nz,
  		p0[1]+nx,p0[2]+ny,p0[3]+nz,
  		0,0, 2,0, 2,0, 0,0,alpha
  		)
  	return 1
  end

  --座標の入ったテーブルtを受け取って線で結ぶ----------------------------------------------------------------------
  -- mode = 0 :順 	maxnum :増えると順に引く、nilで一周。
  -- mode = 1 :すべて巡回 maxnum :増えると順に引く、ネズミ算式に増える。テーブルで第二引数で点あたりの最大数を制限 ・ 省略時は1(mode 0と同じ)
  -- mode = 2 :すべて巡回 maxnum :点一つ当たり何本引くか指定のみ。
  -- tは {{x,y,z},{x,y,z},{x,y,z}...}の形式。 一応 {x,y,z,x,y,z,x,y,z...}でも可。

  local Draw_lineCn = function(mode,t,max_num,w,col,alp,st,va,nst,nva)
  		local cam_param = Camparam()
  		w,col,alp = w or 1, col or 0xffffff, alp or 1
  		if not tostring(col):find("table:") then col = {col} end
          if not tostring(alp):find("table:") then alp = {alp} end
  		st,va,nst,nva = st or 2000, va or 2500,nst or 0,nva or 1
  		if (max_num==nil) then
  			max_num = {#t,1}
  		elseif not tostring(max_num):find("table:") then
  			max_num = {max_num, 1}
  		end
  		local maxnum = max_num[1]
  		local maxcount = max_num[2] or 1
  		if maxcount<0 then maxcount=#t end
  		maxnum,af = math.modf(maxnum+1)
  		af = (math.cos(math.pi*af^2)-1)*-.5
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
  				local color = col[(i-1)%(#col)+1]
                local alpha = alp[(i-1)%(#alp)+1]
  				local progress = 1
  				if maxnum<i then
  					progress = 0
  				elseif maxnum==i then
  					progress = af
  				end
  				n = Draw_line3D(t[i],t[i%#t+1],w,color,alpha,st,va,nil,nil,progress,cam_param)
  			end
        elseif (mode==1) then
            for i=1,#t-1 do
                local n = 0
                for j=i+1,#t do
                    local color = col[(i-1)%(#col)+1]
                    local alpha = alp[(i-1)%(#alp)+1]
                    local progress = 1
                    if maxnum<j then
                        progress = 0
                    elseif maxnum==j then
                        progress = af
                    end
                    --local nc,p1 = Draw_line3D(t[i],t[j%#t+1],w,color,alp,st,va,nil,nil,progress)
                    n = n + Draw_line3D(t[i],t[j],w,color,alpha,st,va,nil,nil,progress,cam_param)
                    if (n>=maxcount) then break end
                end
            end
  		elseif (mode==2) then
            local n = 1
            local tmp = {}
                for i=1,#t-1 do
                    local jn=0
                    for j=i+1,#t do
                        local l = Lensq(t[i],t[j])^.5
                        if (jn>=maxcount) then break end
                        if l>st and l<va then
                            if alp[i]*alp[j]>0 then
                                tmp[n] = {t[i], t[j] ,l=l,col = col[i],alp = alp[i]*alp[j]}
                                jn = jn+1
                                n = n+1
                            end
                        end
                    end
                end
            table.sort(tmp,function(a,b) return (a.l < b.l) end)
    		n = 0
    		for i=1,#tmp do
    				-- local color = col[(i-1)%(#col)+1]
    				-- local alpha = alp[(i-1)%(#alp)+1]
    				local progress = 1
    				if maxnum<i then
    					progress = 0
    				elseif maxnum==i then
    					progress = af
    				end
    				--local nc,p1 = Draw_line3D(t[i],t[j%#t+1],w,color,alp,st,va,nil,nil,progress)
    				n = n + Draw_line3D(tmp[i][1],tmp[i][2],w,tmp[i].col,tmp[i].alp,st,va,nil,nil,progress,cam_param)
    				--if (n>=maxcount) then break end
    		end


  		elseif (mode==3) then
  			-- local maxcount,af = math.modf(maxcount)

  	 		for i=1,#t-1 do
  				local n = 0
  				for j=i,#t do
  					local progress = 1
  					local color = col[n%(#col)+1]
  					local alpha = alp[n%(#alp)+1]
  					--if maxcount<n then
  					--	progress = 0
  					--elseif j==maxcount then
  					--	progress = af
  					--end
  					if Length(Sub(t[i],t[j%#t+1]))<=va then
  						n = n + Draw_line3D(t[i],t[j%#t+1],w,color,alpha,st,va,nil,nil,progress,cam_param)
  					end
  					if (n>=maxnum) then break end
  				end
  			end
  		end
  		return p1
  end
  local Draw_lineCn2 = function(mode,t,max_num,w,col,alp,st,va,nst,nva)
  		local cam_param = Camparam()
  		w,col,alp = w or 1, col or 0xffffff, alp or 1
  		st,va,nst,nva = st or 2000, va or 2500,nst or 0,nva or 1
  		if (max_num==nil) then
  			max_num = {#t,1}
  		elseif not tostring(max_num):find("table:") then
  			max_num = {max_num, 1}
  		end
  		local maxnum = max_num[1]
  		local maxcount = max_num[2] or 1
  		if maxcount<0 then maxcount=#t end
  		maxnum,af = math.modf(maxnum+1)
  		af = (math.cos(math.pi*af^2)-1)*-.5
  -- 		if not tostring(t[1]):find("table:") then
  -- 			local t0={}
  -- 			for i=1,#t/3 do
  -- 				t0[i]={t[i*3-2],t[i*3-1],t[i*3]}
  -- 			end
  -- 			t = t0
  -- 		end
  		--local p1 = t[1] --第二戻り値(線の先端座標)
  		if mode==0 then
  			for i=1,#t do
  				local color = t[i].col
                local alpha = t[i].alpha*alp
  				local progress = 1
  				if maxnum<i then
  					progress = 0
  				elseif maxnum==i then
  					progress = af
  				end
  				n = Draw_line3D(t[i].pos,t[i%#t+1].pos,w,color,alpha,st,va,nil,nil,progress,cam_param)
  			end
        elseif (mode==1) then
            for i=1,#t-1 do
                local n = 0
                for j=i+1,#t do
                    local color = t[i].col
                    local alpha = t[i].alpha*alp
                    local progress = 1
                    if maxnum<j then
                        progress = 0
                    elseif maxnum==j then
                        progress = af
                    end
                    --local nc,p1 = Draw_line3D(t[i],t[j%#t+1],w,color,alp,st,va,nil,nil,progress)
                    n = n + Draw_line3D(t[i].pos,t[j].pos,w,color,alpha,st,va,nil,nil,progress,cam_param)
                    if (n>=maxcount) then break end
                end
            end
  		elseif (mode==2) then
            local n = 1
            local tmp = {}
                for i=1,#t-1 do
                    local jn=0
                    for j=i+1,#t do
                        if (jn>=maxcount) then break end
                        local l = Lensq(t[i].pos,t[j].pos)^.5
                        if l>st and l<(st+va) then
                            if t[i].alpha*t[j].alpha>0 then
                                tmp[n] = {t[i].pos, t[j].pos,l=l,col = t[i].col,alp = t[i].alpha*t[j].alpha}
                                jn = jn+1
                                n = n+1
                            end
                        end
                    end
                end
            --table.sort(tmp,function(a,b) return (a.l < b.l) end)
    		n = 0
    		for i=1,#tmp do
    				local progress = 1
    				if maxnum<i then
    					progress = 0
    				elseif maxnum==i then
    					progress = af
    				end
    				n = n + Draw_line3D(tmp[i][1],tmp[i][2],w,tmp[i].col,tmp[i].alp,st,va,nil,nil,progress,cam_param)
    		end


  		elseif (mode==3) then
  			-- local maxcount,af = math.modf(maxcount)

  	 		for i=1,#t-1 do
  				local n = 0
  				for j=i,#t do
  					local progress = 1
  					local color = t[i].col
  					local alpha = t[i].alpha
  					--if maxcount<n then
  					--	progress = 0
  					--elseif j==maxcount then
  					--	progress = af
  					--end
  					if Length(Sub(t[i].pos,t[j%#t+1].pos))<=va then
  						n = n + Draw_line3D(t[i].pos,t[j%#t+1].pos,w,color,alpha,st,va,nil,nil,progress,cam_param)
  					end
  					if (n>=maxnum) then break end
  				end
  			end
  		end
  		return p1
  end


  -- local Draw_lineCn2 = function(mode,t,max_num,w,col,col2,alp,alp2,st,va,nst,nva)
  -- 		local cam_param = Camparam()
  -- 		w,col,alp = w or 1, col or 0xffffff, alp or 1
  -- 		if not tostring(col):find("table:") then col = {col} end
  --       if not tostring(alp):find("table:") then alp = {alp} end
  -- 		st,va,nst,nva = st or 2000, va or 2500,nst or 0,nva or 1
  -- 		if (max_num==nil) then
  -- 			max_num = {#t,1}
  -- 		elseif not tostring(max_num):find("table:") then
  -- 			max_num = {max_num, 1}
  -- 		end
  -- 		local maxnum = max_num[1]
  -- 		local maxcount = max_num[2] or 1
  -- 		if maxcount<0 then maxcount=#t end
  -- 		maxnum,af = math.modf(maxnum+1)
  -- 		af = (math.cos(math.pi*af^2)-1)*-.5
  -- 		if not tostring(t[1]):find("table:") then
  -- 			local t0={}
  -- 			for i=1,#t/3 do
  -- 				t0[i]={t[i*3-2]-obj.x*0,t[i*3-1]-obj.y*0,t[i*3]-obj.z*0}
  -- 			end
  -- 			t = t0
  -- 		end
  -- 		--local p1 = t[1] --第二戻り値(線の先端座標)
  -- 		if mode==0 then
  -- 			for i=1,#t do
  -- 				local color = col[(i-1)%(#col)+1]
  --               local alpha = alp[(i-1)%(#alp)+1]
  -- 				local color2 = col2[(i-1)%(#col2)+1]
  --               local alpha2 = alp2[(i-1)%(#alp2)+1]
  -- 				local progress = 1
  -- 				if maxnum<i then
  -- 					progress = 0
  -- 				elseif maxnum==i then
  -- 					progress = af
  -- 				end
  -- 				n = Draw_line3D2(t[i],t[i%#t+1],w,color,alpha,st,va,nil,nil,progress,cam_param)
  -- 			end
  --       elseif (mode==1) then
  --           for i=1,#t-1 do
  --               local n = 0
  --               for j=i,#t do
  --                   local color = col[(i-1)%(#col)+1]
  --                   local alpha = alp[(i-1)%(#alp)+1]
  --                   local color2 = col2[(i-1)%(#col2)+1]
  --                   local alpha2 = alp2[(i-1)%(#alp2)+1]
  --                   local progress = 1
  --                   if maxnum<j then
  --                       progress = 0
  --                   elseif maxnum==j then
  --                       progress = af
  --                   end
  --                   --local nc,p1 = Draw_line3D(t[i],t[j%#t+1],w,color,alp,st,va,nil,nil,progress)
  --                   n = n + Draw_line3D2(t[i],t[j],w,color,color2,alpha,alpha2,st,va,nil,nil,progress,cam_param)
  --                   if (n>=maxcount) then break end
  --               end
  --           end
  -- 		elseif (mode==2) then
  --           local n = 1
  --           local tmp = {}
  --               for i=1,#t-1 do
  --                   local jn=0
  --                   for j=i,#t do
  --                       local l = Lensq(t[i],t[j])^.5
  --                       if (jn>=maxcount) then break end
  --                       if l>st and l<va then
  --                           if alp[i]*alp[j]>0 then
  --                               tmp[n] = {t[i], t[j] ,l=l,
  --                                   col=col[i],col2=col2[i],
  --                                   alp=alp[i],alp2=alp2[i]
  --                               }
  --                               jn = jn+1
  --                               n = n+1
  --                           end
  --                       end
  --                   end
  --               end
  --           table.sort(tmp,function(a,b) return (a.l < b.l) end)
  --   		n = 0
  --   		for i=1,#tmp do
  --   				-- local color = col[(i-1)%(#col)+1]
  --   				-- local alpha = alp[(i-1)%(#alp)+1]
  --   				local progress = 1
  --   				if maxnum<i then
  --   					progress = 0
  --   				elseif maxnum==i then
  --   					progress = af
  --   				end
  --   				--local nc,p1 = Draw_line3D(t[i],t[j%#t+1],w,color,alp,st,va,nil,nil,progress)
  --   				n = n + Draw_line3D(tmp[i][1],tmp[i][2],w,tmp[i].col,tmp[i].col2,tmp[i].alp,tmp[i].alp2,st,va,nil,nil,progress,cam_param)
  --   				--if (n>=maxcount) then break end
  --   		end
  --
  --
  -- 		elseif (mode==3) then
  -- 			-- local maxcount,af = math.modf(maxcount)
  --
  -- 	 		for i=1,#t-1 do
  -- 				local n = 0
  -- 				for j=i,#t do
  -- 					local progress = 1
  -- 					local color = col[n%(#col)+1]
  -- 					local alpha = alp[n%(#alp)+1]
  -- 					local color2 = col[n%(#col2)+1]
  -- 					local alpha2 = alp[n%(#alp2)+1]
  -- 					--if maxcount<n then
  -- 					--	progress = 0
  -- 					--elseif j==maxcount then
  -- 					--	progress = af
  -- 					--end
  -- 					if Length(Sub(t[i],t[j%#t+1]))<=va then
  -- 						n = n + Draw_line3D(t[i],t[j%#t+1],w,color,color2,alpha,alpha2,st,va,nil,nil,progress,cam_param)
  -- 					end
  -- 					if (n>=maxnum) then break end
  -- 				end
  -- 			end
  -- 		end
  -- 		return p1
  -- end


--簡易ライティング-------------------------------------------------------------
-- ライトオブジェクトを作成
-- ライトにしたいオブジェクトにスクリプト制御を使いMakelight()
-- 引数は無くても
local Makelight = function(intensity,color,distance)
	if not LightLayer then LightLayer={} end
	intensity	= intensity or 1
	color		= color or obj.getpixel(0,0)
	distance 	= distance or 1200
	distance	= math.max(1,distance)
	local r,g,b 	= RGB(color)
	LightLayer[obj.layer] = {
		intensity = intensity,
		pos = {obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz},
		rot = {-obj.rx,-obj.ry,0},
		distance = distance,
		r = r,
		g = g,
		b = b
		}
	-- アクティブでないライトが他にあったら消しておく(あまり役には立たない)
	local temp = LightLayer
	for k,v in pairs(LightLayer) do
		if k>0 then
			if not obj.getvalue("layer"..k..".x") then
				temp[k]=nil
			end
		end
	end
	LightLayer = temp
end


-- 面重心座標、面法線、カメラの面方向ベクトル
local Face = function(v,camera_param)
	local a={v[1] ,v[2] ,v[3] }
	local b={v[4] ,v[5] ,v[6] }
	local c={v[7] ,v[8] ,v[9] }
	local d={v[10],v[11],v[12]}

	local a_b = Sub(b,a)
	local b_c = Sub(c,b)
	local c_d = Sub(d,c)
	local d_a = Sub(a,d)
	local face0 = Cross(a_b,b_c)
	local face1 = Cross(c_d,d_a)
	local face = Scale( Add(face1, face0) ,.5)
	local surface_normal = Normalize(face)

	-- local maxlen = math.max(Lensq(a,b),Lensq(b,c),Lensq(c,d),Lensq(d,a))
	-- maxlen = math.sqrt(maxlen)
	local center = {(a[1]+b[1]+c[1]+d[1])/4,(a[2]+b[2]+c[2]+d[2])/4,(a[3]+b[3]+c[3]+d[3])/4}
    local cam = obj.getoption("camera_param")
	local eye = Normalize(  Sub(center,{cam.x,cam.y,cam.z} ))
	return center,surface_normal,eye
end

--require("rikky_module")
local Refx = function (vertex, intensity, material_col, ambient_col, lightlayer , blend, min)
	-- ライトテーブルが見つからなければ何もしない
	if not lightlayer or #lightlayer<1 then
        return
    end

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
	material_col = material_col or obj.getpixel(0,0)
	local R,G,B    = RGB(material_col)
	obj.pixeloption("type","rgb")
	obj.putpixel(0,0,R,G,B,255)
	R,G,B = R/255,G/255,B/255
	--local ambient_col = light_col[2]
	local aR,aG,aB = RGB(ambient_col)
	local aR,aG,aB = aR*R*intensity_ambient, aG*G*intensity_ambient, aB*B*intensity_ambient
	-- // rikky_moduleを使う場合 //
	-- local cam = rikky_module.camerainfo(vertex)		 -- vertexはdrawpoly用の頂点4セット
	-- local surface = {cam.mx, cam.my, cam.mz}		 -- 面法線
	-- local eye = {cam.vx, cam.vy, cam.vz}			 -- 面へのカメラ視線ベクトル
	-- local center = Face(vertex)				 -- 面の重心座標

	-- // rikky_moduleがない場合 //
	local center,surface,eye = Face(vertex)
    center = {Rot(center)}
	center = Add(center,{obj.x, obj.y, obj.z})		 -- 基準座標のみ追加
	local vc = Dot(eye ,surface)				 -- 面がカメラ方向に向いているか (カメラ==ライトの場合はこれだけ)
	local dR,dG,dB = 0,0,0
	local rR,rG,rB = 0,0,0
	local mR,mG,mB = R,G,B
	--obj.putpixel(0,0,R,G,B,255)
	local rt = {}
	for k,light in pairs(lightlayer) do
		local lR,lG,lB = light.r, light.g, light.b
		lR,lG,lB = lR/255,lG/255,lB/255
		local lcol = {lR,lG,lB}
		local light_vec,dist = Normalize(Sub(light.pos,center))	 -- 面→ライトのベクトル、面との距離
        light_vec = {Rot(light_vec,light.rot)} 						--ライトの回転
		local vd0 = Dot(surface,light_vec)	 		 -- vd,面がライトに向いているか(単純なDiffuse用)
		local inter_eye = (Dot(eye,surface) * vd0)>0 and 0 or 1	 -- 光源とカメラがポリゴンで遮られているか簡易判定
		vd0 = math.max(0,math.min(1,vd0)) --* inter_eye
		local power = 1/math.max(1,(dist/light.distance)^2)		 -- 距離減衰
		local vd = vd0 * intensity_diffuse * power * light.intensity
		local vr = Vec_reflect(surface ,light_vec)		 -- vr,面に反射した後のライトのベクトル
		local v2 = Dot(eye ,vr)				 	 -- v2,反射光がカメラに向いているか(反射光用)
		v2 = math.max(0,v2) * inter_eye
		v2 = Linear(v2,min,1,0,1) * intensity_specular * power
		V2_RET = v2

		dR,dG,dB = dR*1 + aR*0 + 255*mR*lR*vd,  dG*1 + aG*0 + 255*mG*lG*vd,  dB*1 + aB*0 + 255*mB*lB*vd
		rR,rG,rB = (dR*1+rR + 255*lR)*v2, (dG*1+rG + 255*lG)*v2, (dB*1+rB + 255*lB)*v2
		rR,rG,rB = math.min(rR,255), math.min(rG,255), math.min(rB,255)
		rt[k] = {r=rR, g=rG ,b=rB, a=v2*255* light.intensity}				 -- スペキュラはライト個別で加算
	end
	obj.pixeloption("blend",0)
	obj.putpixel(0,0,math.min(255,dR+aR),math.min(255,dG+aG),math.min(255,dB+aB),255)

	-- スペキュラ
	for k,v in pairs(rt) do
		obj.pixeloption("blend",blend or 0)
		obj.putpixel(0,0, v.r, v.g, v.b, math.min(255,v.a))
	end
	obj.pixeloption("blend")
	obj.pixeloption("type","col")
    return obj.getpixel(0,0)
end

local Maxlen = function(a,b,c,d)
    local v
    if d then
        v=math.max(Lensq(a,b),Lensq(b,c),Lensq(c,d),Lensq(d,a))
    else
        v=math.max(Lensq(a,b),Lensq(b,c),Lensq(c,a))
    end
    return math.sqrt(v)
end

--本体 ----------------------------------------------------------------------------------
local FacetsXY = function(t,col,alp,tri,reflect,start_vanish,vanish)
    local param = Param()
    -- obj.ry,obj.rx,obj.rz=0,0,0
	local x,y,z = 0,0,0 --Rotc() --obj.x+obj.ox, obj.y+obj.oy, obj.z+obj.oz
	obj.setoption("billboard",0)
	alp = alp or 1
	tri = tri or 0
    -- local alpha = alp
    start_vanish,vanish = start_vanish or 500,vanish or 3000
    vanish = start_vanish+vanish
	local Ref = function() return end
	local Put = function() return end
	reflect = reflect or {}
	if #reflect>0 then
		Ref = Refx
	else
		reflect = {0,nil,nil,nil,nil,0}
	end
	obj.setoption("billboard",0)
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
	tri = tri or 0
	--if col then obj.putpixel(0,0,col,1) end
	local xn=#t
	local yn=#t[1]
	local zn=#t[1][1]
	if (tri==0) then
		for i=1,xn-1 do
		for j=1,#t[1]-1 do
		for k=1,#t[i][j] do
            local t0=t[i][j][k]
			local t1=t[i+1][j][k]
			local t2=t[i+1][j+1][k]
			local t3=t[i][j+1][k]
			local x0,y0,z0 = t0[1]	,t0[2]	,t0[3]
			local x1,y1,z1 = t1[1]	,t1[2]	,t1[3]
			local x2,y2,z2 = t2[1]	,t2[2]	,t2[3]
			local x3,y3,z3 = t3[1]	,t3[2]	,t3[3]
            local al,Y = 1,1
            if #t[i][j][1]==8 then
                local C = t[i][j][1]
                local Ca = t[i+1][j+1][1]
                local r,g,b = C[4],C[5],C[6]
                al,Y = C[7]*Ca[7],C[8]
                local COL = RGB(r,g,b)
                reflect[2] = COL
                col = {COL}
            end
            local alp=alp*al
            local maxlen = Maxlen(t0,t1,t2,t3)
            local alpha = alp * Linear(maxlen,start_vanish,vanish,1,0)
            if (alpha>0) then
                Put(0,0,col[i%(#col)+1],1)
                Ref({x0+x,y0+y,z0+z, x1+x,y1+y,z1+z, x2+x,y2+y,z2+z, x3+x,y3+y,z3+z},unpack(reflect))
                obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha)
            end
		end
		end
		end
	elseif (tri==1) then
		for i=1,xn-1 do
		for j=1,#t[1]-1 do
		for k=1,#t[i][j] do
            local t0=t[i][j][k]
			local t1=t[i+1][j][k]
			local t2=t[i+1][j+1][k]
			local t3=t[i][j+1][k]
			local x0,y0,z0 = t0[1]	,t0[2]	,t0[3]
			local x1,y1,z1 = t1[1]	,t1[2]	,t1[3]
			local x2,y2,z2 = t2[1]	,t2[2]	,t2[3]
			local x3,y3,z3 = t3[1]	,t3[2]	,t3[3]
            local al,Y = 1,1
            if #t[i][j][1]==8 then
                local C0 = t[i][j][1]
                local C1 = t[i+1][j][1]
                local C2 = t[i+1][j+1][1]
                local C3 = t[i][j+1][1]
                local r,g,b = C0[4],C0[5],C0[6]
                al,Y = C0[7]*C1[7]*C3[7]*C3[7], C0[8]
                local COL = RGB(r,g,b)
                reflect[2] = COL
                col = {COL}
            end
            local alp=alp*al
			if (j%2==0) then
                local maxlen = Maxlen(t0,t1,t2)
                local alpha = alp * Linear(maxlen,start_vanish,vanish,1,0)
                if (alpha>0) then
                    Put(0,0,col[i%(#col)+1],1)
                    Ref({x0+x,y0+y,z0+z, x1+x,y1+y,z1+z, x2+x,y2+y,z2+z, x2+x,y2+y,z2+z},unpack(reflect))
                    obj.drawpoly(x0,y0,z0, x1,y1,z1, x2,y2,z2, x2,y2,z2, 0,0,0,0,0,0,0,0,alpha)
                end
                 maxlen = Maxlen(t0,t2,t3)
                 alpha = alp * Linear(maxlen,start_vanish,vanish,1,0)
                 if (alpha>0) then
                     Put(0,0,col[i%(#col)+1],1)
                     Ref({x0+x,y0+y,z0+z, x2+x,y2+y,z2+z, x3+x,y3+y,z3+z, x3+x,y3+y,z3+z},unpack(reflect))
                     obj.drawpoly(x0,y0,z0, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha)
                 end
			else
                local maxlen = Maxlen(t0,t1,t3)
                local alpha = alp * Linear(maxlen,start_vanish,vanish,1,0)
                if (alpha>0) then
                    Put(0,0,col[i%(#col)+1],1)
                    Ref({x0+x,y0+y,z0+z, x1+x,y1+y,z1+z, x3+x,y3+y,z3+z, x3+x,y3+y,z3+z},unpack(reflect))
                    obj.drawpoly(x0,y0,z0, x1,y1,z1, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha)
                end
                maxlen = Maxlen(t1,t2,t3)
                alpha = alp * Linear(maxlen,start_vanish,vanish,1,0)
                if (alpha>0) then
                    Put(0,0,col[i%(#col)+1],1)
                    Ref({x1+x,y1+y,z1+z, x2+x,y2+y,z2+z, x3+x,y3+y,z3+z, x3+x,y3+y,z3+z},unpack(reflect))
                    obj.drawpoly(x1,y1,z1, x2,y2,z2, x3,y3,z3, x3,y3,z3, 0,0,0,0,0,0,0,0,alpha)
                end
			end
		end
		end
		end
	end

end

local Depthfx = function(
	pos,		-- オブジェクト座標{x,y,z}
	focalpoint,	-- 焦点の前後
	startfade,	-- フェード開始位置
	vanish,		-- フェード開始からの消滅距離
	near_startfade, -- フェード開始位置(焦点より手前)
	near_vanish,	-- フェード開始からの消滅距離 (焦点より手前)
	focusmode,	-- trueで焦点を目標点に固定
    blur,		-- ﾌﾞﾗｰ強度
    aspect,     --
	blur_mode,	--0 or nil = ぼかし、1=レンズブラー
	alpha,		-- 1で透明度変化
	cam_param,  -- Camparamの戻り値をあらかじめ入れる用
	depth_map  -- 既に使用しているDepthfxの戻り値、エフェクトのみ、計算しない
	)
	pos = pos or {obj.x+obj.ox,obj.y+obj.oy,obj.z+obj.oz}
	focalpoint,startfade,vanish = focalpoint or 0,startfade or 100,vanish or 2500
	near_startfade,near_vanish = near_startfade or startfade,near_vanish or vanish
	alpha = alpha or 0
    aspect = aspect or 0
    local nvanish = 100
	local D,depth,pd,nD=0,0,0,1
	if not depth_map then
		local cam = cam_param --or obj.getoption("camera_param")
		local e = cam.eye --or {c.tx-c.x, c.ty-c.y, c.tz-c.z}
		local d = cam.l --or (e[1]*e[1]+e[2]*e[2]+e[3]*e[3])^.5
		local n = cam.n --or Scale(e,1/l)
		local pd = (focusmode) and (d+focalpoint) or (1024+focalpoint)
		local pl  = Scale(n,pd)
		local pv = {pos[1]-cam.x, pos[2]-cam.y, pos[3]-cam.z}
		local depth  = Dot({pv[1]-pl[1],pv[2]-pl[2],pv[3]-pl[3]},n)
        local npl = Scale(n,nvanish)
        local ndepth = Dot({pv[1]-npl[1],pv[2]-npl[2],pv[3]-npl[3]},n)
        if ndepth<0 then
            nD = math.abs(ndepth)-nvanish
            nD = nD<0 and 0 or nD
            nD = nD>nvanish and nvanish or nD
            nD = (1-nD/nvanish)
        end
		if depth<0 then
            startfade,vanish = near_startfade,near_vanish
        end
		D = math.abs(depth)-startfade
		D = D<0 and 0 or D
		D = D>vanish and vanish or D
		D = (1-D/vanish)

	else
		D,depth,pd,nD = depth_map[1],depth_map[2],depth_map[3],depth_map[4]
	end
	if alpha==1 then obj.alpha = obj.alpha*D end
	if (obj.alpha>0) then
	    if blur>0 then
	        if (blur_mode==1) then
	            obj.effect("レンズブラー","サイズ固定",0,"範囲",blur*(1-D))
	        else
	            obj.effect("ぼかし","サイズ固定",0,"範囲",blur*(1-D),"縦横比",aspect)
	        end
	    end
	end
	return {D,depth,pd,nD}
end

Vector =  {
    Sub    = Sub,
    Add    = Add,
    Scale  = Scale,
    Len    = Len,
    Lensq  = Lensq,
    Length = Length,
    Norm   = Normalize,
    Normalize = Normalize,
    Cross  = Cross,
    Dot    = Dot,
    Mul    = Mul,
    Div    = Div,

    Dist_p_line    = Dist_p_line,
    Dist_p_segment = Dist_p_segment,
    Pos_p_on_plane = Pos_p_on_plane,
    Pos_p_on_poly  = Pos_p_on_poly,
    Pos_p_on_line                  = Pos_p_on_line,
    Pos_p_on_segment               = Pos_p_on_segment,
    Pos_plane_intersection_segment = Pos_plane_intersection_segment,
    Pos_intersection_2line         = Pos_intersection_2line,

    Vec_reflect = Vec_reflect,
    Norm_cameraplane = Norm_cameraplane,
    Norm_surface = Norm_surface
}

Drawline = {
      P2 = Draw_line,
      P3 = Draw_line3D,
      Pc = Draw_lineCn
  }

getcolortools ={
    Vector = Vector,
    Linear = Linear,
    TBL = TBL,
    GL = GL,
    GL2= GL2,
    Param = Param,
    Camparam = Camparam,
    Progress = Progress,
    Shuffle = Shuffle,
    Linear = Linear,
    Swap_table = Swap_table,
    Swap_i_table = Swap_i_table,
    Swap_col = Swap_col,
    Rot = Rot,
    Rot_c = Rot_c,
    Rot_rpy = Rot_rpy,
    Rot_v = Rot_v,
    Co_spherical = Co_spherical,
    Co_torus = Co_torus,
    Makelight = Makelight,
    FacetsXY = FacetsXY,
    Draw_line = Draw_line,
    Draw_line3D = Draw_line3D,
    Draw_line3D2 = Draw_line3D2,
    Draw_lineCn = Draw_lineCn,
    Draw_lineCn2 = Draw_lineCn2,
    Depthfx = Depthfx
}
return false
